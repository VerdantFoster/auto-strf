#!/usr/bin/env python3
"""Compatibility entrypoint: uses the maintained calibrated waterfall helper."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from satnogs_waterfall.__main__ import main
if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1].isdigit():sys.argv.insert(1,'mark')
    main()
