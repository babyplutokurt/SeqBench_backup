import sys

sys.path.append('./')
sys.path.append('./Compression_Scripts')
sys.path.append('./Post_Hoc_Scripts')
from Compression_Scripts.job_generator import JobGenerator
from Post_Hoc_Scripts.post_hoc_analysis import PostHocAnalysis
from Error_Analysis_Scripts.error_analysis import ErrorAnalysis

if __name__ == '__main__':
    config = '/home/tus53997/SeqBench/Jobs/bench_HG0097.json'
    jb = JobGenerator(config)
    jb.generate_and_submit_jobs()
    ea = ErrorAnalysis(config)
    ea.run_error_analysis()
    pa = PostHocAnalysis(config)
    pa.run_posthoc_analysis()
