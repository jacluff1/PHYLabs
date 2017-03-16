# Run noise.py
import noise2 as no
from importlib import reload
reload(no)

parts = [ no.part2() , no.part3() ]

runAll=False

if runAll == True:
    for part in parts:
        part
else:
    parts[0]
