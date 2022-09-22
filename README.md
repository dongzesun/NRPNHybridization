# HybridizationWaveforms

The current version do hybridization between CCE waveforms and PN waveforms.
```
python Hybridization.py --t --SimDir --CCEDir --OutDir --length
    '--t',type=float, default=-7000.0,help='Start time of matching window'
    '--SimDir', default='/panfs/ds09/sxs/dzsun/SimAnnex/Public/HybTest/015/Lev3',help='Path in which to find the extropolated waveform data'
    '--CCEDir', default='/home/dzsun/CCEAnnex/Public/HybTest/015_CCE/Lev3/CCE',help='Path in which to find the CCE waveform data'
    '--OutDir', default='/home/dzsun',help='Path in which to output results'
    '--length',type=float, default=5000.0,help='Length of matching window'
    '--nOrbits',type=float, default=None,help='Length of matching window in orbits, will disable "length" option if not None'
    '--truncate',nargs=2,type=float, default=None,help='--truncate t1 t2. If specified, it will truncate the abd object and keep only data between t1 and t2'
```
Note: if the datasize is too large, please enable '--truncate' option by adding `--truncate t1 t2` argument, so that the data is truncated and will use less memory.
