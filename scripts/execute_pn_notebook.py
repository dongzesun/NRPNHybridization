#!/usr/bin/env python3
"""Execute legacy symbolic notebooks with explicit, fail-fast Python 3 semantics.

The historical ExecNotebook.ipy depends on Python 2 execfile and deprecated
IPython notebook APIs, and ignores errors in nested cells. Supply its two
execution helpers directly, keeping the source notebooks and shared namespace.
"""
import ast
import json
import os
from pathlib import Path
import sys


def run(path):
    namespace={'__name__':'__main__','display':lambda *args,**kwargs:None}

    def code(source,filename):
        tree=ast.parse(source)
        # These notebooks place the obsolete division future after imports.
        # Python 3 already supplies the intended true-division semantics.
        tree.body=[n for n in tree.body if not (isinstance(n,ast.ImportFrom) and n.module=='__future__')]
        exec(compile(tree,filename,'exec'),namespace)

    def execfile(filename):
        if Path(filename).name=='ExecNotebook.ipy':
            return
        code(Path(filename).read_text(),str(Path(filename).resolve()))

    def execnotebook(filename):
        filename=Path(filename).resolve()
        previous=Path.cwd()
        os.chdir(filename.parent)
        try:
            nb=json.loads(filename.read_text())
            cells=nb['cells'] if 'cells' in nb else nb['worksheets'][0]['cells']
            for cell in cells:
                if cell['cell_type']=='code':
                    code(''.join(cell.get('source',cell.get('input',[]))),str(filename))
        finally:
            os.chdir(previous)

    namespace.update(execfile=execfile,execnotebook=execnotebook)
    execnotebook(path)


if __name__=='__main__':
    run(Path(sys.argv[1]))
