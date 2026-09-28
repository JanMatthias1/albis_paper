"""Accounting checks: endpoint semantics, shared reference cost, and failed runs."""
import json
from pathlib import Path
import tempfile
import unittest
from native_worker import Measurements, save
from summarize_native import collect, review_pilot


class NativeAccountingTests(unittest.TestCase):
    def test_reference_once_and_validation_excluded(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            save(root/'protocol.json',dict(protocol_id='test',smoke_test=True))
            task=dict(key='spider_n32_seed1',method='spider',n_cells=32,seed=1)
            save(root/'tasks.json',[task])
            out=root/'raw'/task['key'];(out/'worker').mkdir(parents=True)
            save(out/'process_measurement.json',dict(returncode=0,wall_seconds=100,peak_tree_rss_bytes=1000))
            stage=lambda name,category,seconds,rss:dict(name=name,category=category,status='ok',wall_seconds=seconds,cpu_seconds=seconds/2,peak_tree_rss_bytes=rss)
            measurement=dict(status='ok',endpoint_validated=True,stages=[stage('imports','startup',1,100),stage('load','input',2,200),stage('spatial','generation',3,300),stage('expression','generation',4,400),stage('validate','validation',90,900)],
                reference_measurement=dict(wall_seconds=5,peak_tree_rss_bytes=500),output=dict(n_cells=32,n_genes=556))
            save(out/'worker/measurement.json',measurement)
            _,rows,_=collect(root,False);row=rows[0]
            self.assertEqual(row['generation_seconds'],7)
            self.assertEqual(row['native_pipeline_seconds'],10)
            self.assertEqual(row['from_scratch_seconds'],15)
            self.assertEqual(row['from_scratch_peak_gib'],500/2**30)
            self.assertEqual(row['generation_peak_gib'],400/2**30)
            measurement['status']='failed';save(out/'worker/measurement.json',measurement)
            _,rows,_=collect(root,False)
            self.assertIsNone(rows[0]['generation_seconds'])
            self.assertEqual(rows[0]['status'],'failed')

    def test_stage_failure_is_persisted(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp);report=dict(stages=[])
            measure=Measurements(out,report)
            try:
                with self.assertRaisesRegex(ValueError,'sentinel'):
                    with measure.stage('native','generation'):
                        raise ValueError('sentinel')
            finally: measure.close()
            saved=json.loads((out/'measurement.json').read_text())['stages'][0]
            self.assertEqual(saved['status'],'failed')
            self.assertGreater(saved['peak_tree_rss_bytes'],0)
            self.assertGreaterEqual(saved['rss_samples'],2)

    def test_pilot_failure_blocks_full_sweep(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            review_pilot(root,dict(smoke_test=False),[])
            self.assertEqual(json.loads((root/'pilot_review.json').read_text())['status'],'incomplete')

class PilotReviewTests(unittest.TestCase):
    def test_complete_pilot_resources_and_hardware_gate(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            cfg=dict(smoke_test=False,sizes=[10000,100000,1000000])
            rows=[dict(method=method,seed=2025,n_cells=n,status='ok',cpu_model='same CPU',
                       process_wall_seconds=n/100,process_peak_gib=n/10000)
                  for method in ['albis','sccube','spider'] for n in [10000,100000]]
            review_pilot(root,cfg,rows)
            report=json.loads((root/'pilot_review.json').read_text())
            self.assertEqual(report['status'],'passed')
            self.assertEqual(report['resources']['albis']['1000000']['memory'],'204G')
            rows[0]['cpu_model']='different CPU'
            review_pilot(root,cfg,rows)
            self.assertEqual(json.loads((root/'pilot_review.json').read_text())['status'],'needs_review')

if __name__=='__main__': unittest.main()
