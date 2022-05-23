from pathlib import Path
from typing import Annotated, List, Optional, Tuple
import sys
sys.path.insert(0, '/root/wf/')
import SpacerMakerCas13b
from latch import small_task, workflow
from latch.types import LatchDir, LatchFile

@small_task
def SpacerMakerCas13b_tsk(genbank_file: LatchFile, run_name: str) -> (LatchFile):
    output_file = SpacerMakerCas13b.execute(genbank_file.local_path, run_name)
    return LatchFile(output_file, "latch:///{}".format(run_name))

@workflow
def SpacerMakerCas13b_wf(genbank_file: LatchFile, run_name: str) -> (LatchFile):
    """
    __metadata__:
        display_name: SpacerMakerCas13b
        author: Matthew Johnson
            name: SpacerMakerCas13b
            email: matthew.johnson@ucsf.edu
            github: https://github.com/mjoh223
        repository: https://github.com/mjoh223/SpacerMakerCas13b
        license:
            id: MIT
    Args:
        genbank_file:
          input gb file.
          __metadata__:
            display_name: genbank file
        run_name:
          run name
          __metadata__:
            display_name: run name
    """
    return SpacerMakerCas13b_tsk(genbank_file = genbank_file, run_name = run_name)
