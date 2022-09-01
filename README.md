# HybridizationWaveforms

The current version do hybridization between CCE waveforms and PN waveforms.
```
python Hybridization.py --t --SimDir --CCEDir --OutDir --length
    '--t',type=float, default=-7000.0,help='Start time of matching window'
    '--SimDir', default='/panfs/ds09/sxs/dzsun/SimAnnex/Public/HybTest/015/Lev3',help='Path in which to find the extropolated waveform data'
    '--CCEDir', default='/home/dzsun/CCEAnnex/Public/HybTest/015_CCE/Lev3/CCE',help='Path in which to find the CCE waveform data'
    '--OutDir', default='/home/dzsun',help='Path in which to output results'
    '--length',type=float, default=5000.0,help='Length of matching window'
```
