#!/usr/bin/env python3

import sys
from unittest.mock import MagicMock

mock_modules = {
    "SigProfilerAssignment": MagicMock(),
    "SigProfilerAssignment.Analyzer": MagicMock(),
    "SigProfilerMatrixGenerator": MagicMock(),
    "SigProfilerMatrixGenerator.scripts": MagicMock(),
}

for name, mock in mock_modules.items():
    sys.modules[name] = mock
