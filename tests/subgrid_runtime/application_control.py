"""Staged matched application test; run only inside a separately admitted build.

Arguments: solver, byte-exact public native XML, fresh output directory,
cpu|reference|reference-two|cuda|cuda-two. No launch/admission or device selection is implicit.
"""
from pathlib import Path
import hashlib,importlib.util,json,os,re,subprocess,sys,time,xml.etree.ElementTree as ET
import numpy as np
import h5py

SOLVER=Path(sys.argv[1]).resolve();XML=Path(sys.argv[2]).resolve();OUT=Path(sys.argv[3]).resolve();MODE=sys.argv[4]
assert MODE in ['cpu','reference','reference-two','cuda','cuda-two','ordinary-cpu','ordinary-reference']
assert 'dipole' in XML.read_text() and 'CartesianSubgrid' in XML.read_text()
OUT.mkdir(exist_ok=False)
input_bytes=XML.read_bytes()
if MODE.startswith('ordinary-'):
 tree=ET.fromstring(input_bytes);fdtd=tree.find('FDTD');assert fdtd is not None
 keys=[k for k in fdtd.attrib if k.startswith('CartesianSubgrid')]
 assert set(keys)=={'CartesianSubgridBox','CartesianSubgridRatio'}
 for key in keys:fdtd.attrib.pop(key)
 input_bytes=ET.tostring(tree,encoding='utf-8',xml_declaration=True)

SOURCE=Path(os.environ['SUBGRID_TEST_SOURCE'])
def load(name):
 spec=importlib.util.spec_from_file_location(name,SOURCE/'python/openEMS'/(name+'.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
control=load('run_control');snapshot=load('snapshots');startup=load('startup_gate')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
flags=['--engine=basic' if MODE in ['cpu','ordinary-cpu'] else '--engine=gpu-reference' if (MODE.startswith('reference') or MODE=='ordinary-reference') else '--engine=gpu','--numThreads=1','--exact-endcriteria']
workers=0 if MODE.startswith('ordinary-') else 1
if MODE in ['cuda-two','reference-two']:
 flags+=['--subgrid-domain-axis=x','--subgrid-domain-cuts=0','--subgrid-domain-workers='+('reference,reference' if MODE=='reference-two' else 'auto@0,auto@1')];workers=2
receipt=dict(status='RUNNING',solverSHA256=sha(SOLVER),inputSHA256=hashlib.sha256(input_bytes).hexdigest(),mode=MODE,controls=[])
def save():(OUT/'control-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
def wait(predicate,child,description,timeout=120):
 deadline=time.monotonic()+timeout
 while time.monotonic()<deadline:
  result=predicate()
  if result:return result
  assert child.poll() is None,(description,child.returncode)
  time.sleep(.02)
 raise TimeoutError(description)
def witness(path,tag,token):
 if not path.exists():return None
 root=ET.parse(path).getroot()
 if root.get('state')!='paused' or root.get('id')!=token:return None
 assert root.tag==tag and int(root.get('nativeCompletedWorkers'))==workers
 assert root.get('nativeCompletedEpoch')==root.get('timestep')
 return dict(root.attrib)
def read_hdf(path):
 data={}
 for f in sorted(path.glob('*.h5')):
  with h5py.File(f)as h:
   def visit(n,v):
    if isinstance(v,h5py.Dataset):data[f.name+'/'+n]=v[...]
   h.visititems(visit)
 return data
def read_arrays(path):
 data=read_hdf(path)
 for n in ['port_ut_1','port_it_1','arm_lower_i','arm_upper_i']:
  data[n]=np.loadtxt(path/n,comments='%')
 assert data
 return data
try:
 for controlled in [False,True]:
  p=OUT/('paused' if controlled else 'control');p.mkdir();(p/'input.xml').write_bytes(input_bytes)
  args=[str(SOLVER),'input.xml',*flags]
  if controlled:args+=['--startup-gate=132','--startup-gate-id=probe-start','--startup-gate-timeout=120','--runtime-control','--decay-gate-id=probe-run','--decay-gate-timeout=120','--dft-snapshots','--decay-gates=-10,-20']
  with(p/'solver.log').open('x')as log:
   child=subprocess.Popen(args,cwd=p,stdout=log,stderr=subprocess.STDOUT)
   try:
    if controlled:
     first=wait(lambda:witness(p/'startup-gate.xml','StartupGate','probe-start'),child,'startup acknowledgement')
     assert int(first['timestep'])>=132
     target=snapshot.request_dft_snapshot(p,'s-start')
     wait(lambda:(target/'snapshot.xml').is_file(),child,'snapshot while startup-paused')
     time.sleep(.1);assert witness(p/'startup-gate.xml','StartupGate','probe-start')['timestep']==first['timestep']
     repeated=snapshot.request_dft_snapshot(p,'s-start-repeat')
     wait(lambda:(repeated/'snapshot.xml').is_file(),child,'repeat snapshot at same paused epoch')
     first_arrays=read_hdf(target);again=read_hdf(repeated);assert first_arrays.keys()==again.keys() and first_arrays
     for name in first_arrays:np.testing.assert_array_equal(first_arrays[name],again[name])
     for f in target.glob('*.h5'):
      with h5py.File(f)as x,h5py.File(repeated/f.name)as y:
       for key in ['sample_count','last_sample_timestep','last_sample_time','sample_interval']:
        assert x.attrs[key]==y.attrs[key],(f.name,key)

     # Queue the next pause while startup remains stopped: no timing race with a
     # fast solver completing before a second request can be submitted.
     control.control_run(p,'pause','probe-run-manual')
     startup.decide_startup_gate(p,'probe-start',True)
     second=wait(lambda:witness(p/'runtime-gate.xml','RuntimeGate','probe-run-manual'),child,'manual acknowledgement')
     assert int(second['timestep'])>=int(first['timestep'])
     target2=snapshot.request_dft_snapshot(p,'s-manual')
     wait(lambda:(target2/'snapshot.xml').is_file(),child,'snapshot while manual-paused')
     time.sleep(.1);assert witness(p/'runtime-gate.xml','RuntimeGate','probe-run-manual')['timestep']==second['timestep']
     for root,ack in [(target,first),(target2,second)]:
      m=ET.parse(root/'snapshot.xml').getroot();assert m.get('nativeCompletedEpoch')==ack['timestep'] and int(m.get('nativeCompletedWorkers'))==workers
      for f in root.glob('*.h5'):
       with h5py.File(f)as h:
        assert h.attrs['sample_count']>0 and int(h.attrs['last_sample_timestep'])<=int(ack['timestep'])
        assert not bool(h.attrs['simulation_converged'])
     for f in target.glob('*.h5'):
      with h5py.File(f)as x,h5py.File(target2/f.name)as y:
       assert y.attrs['sample_count']>=x.attrs['sample_count']
       assert y.attrs['last_sample_timestep']>=x.attrs['last_sample_timestep']
     control.control_run(p,'continue','probe-run-manual')
     # A queued manual request may be acknowledged immediately after startup,
     # at the same complete epoch. Later decay gates establish resumed progress.
     acknowledgements=[first,second]
     signal=re.search(r'Excitation signal length is: (\d+)',(p/'solver.log').read_text())
     assert signal is not None
     for level in [10,20]:
      token='probe-run-d'+str(level)
      ack=wait(lambda:witness(p/'runtime-gate.xml','RuntimeGate',token),child,'post-pulse decay acknowledgement')
      assert ack['reason']=='decay' and int(ack['thresholdDb'])==-level
      assert float(ack['energyDb'])<=-level and int(ack['timestep'])>=int(signal[1])
      assert int(ack['timestep'])>=int(acknowledgements[-1]['timestep'])
      root=p/'dft-snapshots'/token
      assert (root/'snapshot.xml').exists() and ack['snapshotAvailable']=='true'
      for f in root.glob('*.h5'):
       with h5py.File(f)as h:
        assert h.attrs['sample_count']>0 and int(h.attrs['last_sample_timestep'])<=int(ack['timestep'])
      acknowledgements.append(ack)
      control.control_run(p,'continue',token)
     receipt['controls']=acknowledgements;save()
    assert child.wait(timeout=300)==0,(p,child.returncode)
   finally:
    if child.poll()is None:
     child.terminate()
     try:child.wait(timeout=10)
     except subprocess.TimeoutExpired:child.kill();child.wait()
  assert sha(p/'input.xml')==receipt['inputSHA256']
 a=read_arrays(OUT/'control');b=read_arrays(OUT/'paused');assert a.keys()==b.keys()
 for n in a:np.testing.assert_array_equal(a[n],b[n])
 # Native parsing must also reject a forged controller file rather than
 # silently resuming a mismatched decision. These are new fresh subprocesses.
 for fault in ['manual-id','startup-decision']:
  p=OUT/('negative-'+fault);p.mkdir();(p/'input.xml').write_bytes(input_bytes)
  args=[str(SOLVER),'input.xml',*flags,'--startup-gate=132','--startup-gate-id=probe-start','--startup-gate-timeout=120','--runtime-control','--decay-gate-id=probe-run','--dft-snapshots']
  with(p/'solver.log').open('x')as log:
   child=subprocess.Popen(args,cwd=p,stdout=log,stderr=subprocess.STDOUT)
   try:
    wait(lambda:witness(p/'startup-gate.xml','StartupGate','probe-start'),child,'negative startup acknowledgement')
    if fault=='manual-id':
     (p/'RUNTIME_PAUSE').write_text('foreign-token\n')
     startup.decide_startup_gate(p,'probe-start',True)
    else:
     (p/'STARTUP_CONTINUE').write_text('foreign-token\n')
    assert child.wait(timeout=30)==1,(fault,child.returncode)
   finally:
    if child.poll()is None:
     child.terminate()
     try:child.wait(timeout=10)
     except subprocess.TimeoutExpired:child.kill();child.wait()
  raw=(p/'solver.log').read_text()
  assert ('invalid manual pause ID' if fault=='manual-id' else 'decision ID mismatch')in raw
  assert not(p/'runtime-gate.xml').exists(),'invalid request published a valid runtime pause'
 receipt['nativeControllerNegatives']=['manual-id','startup-decision']
 assert sha(SOLVER)==receipt['solverSHA256']
 receipt.update(status='PASS_NATIVE_APPLICATION_PAUSE_SNAPSHOT_RESUME',exactFinalArrays=len(a),workers=workers,limits=['Small public dipole only','No Metal/physical acceptance unless actual admitted mode confirms it','Timeout and injected backend-failure application controls remain separate; core worker-failure gates are R60'])
except BaseException as e:receipt.update(status='FAILED_NO_RESTART',error=repr(e));raise
finally:save()
print(receipt['status'])
