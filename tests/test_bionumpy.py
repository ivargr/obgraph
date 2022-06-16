import pytest
import time
import numpy as np
from bionumpy.delimited_buffers import VCFBuffer, VCFMatrixBuffer, PhasedVCFMatrixBuffer
import bionumpy as bnp


@pytest.mark.skip
def test():
    t = time.perf_counter()
    file = bnp.open("variants.vcf", chunk_size=10000000, buffer_type=PhasedVCFMatrixBuffer, mode="stream")
    n = 0
    for chunk in file:
        print(np.unique(chunk.genotypes))
        n += chunk.genotypes.shape[0]
        print(n)
        print(time.perf_counter()-t)

test()