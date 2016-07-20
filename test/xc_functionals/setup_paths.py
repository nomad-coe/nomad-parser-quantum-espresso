import sys, os, os.path
baseDir = os.path.dirname(os.path.abspath(__file__))
commonDir = os.path.normpath(os.path.join(baseDir,"../../../../python-common/common/python"))
qeDir = os.path.normpath(os.path.join(baseDir,"../../parser/parser-quantum-espresso"))

if not commonDir in sys.path:
    sys.path.insert(0, commonDir)
if not qeDir in sys.path:
    sys.path.insert(0, qeDir)
