import subprocess
import sys
import os

srcdir = "/".join(__file__.split('/')[:-1]) + "/src"
print(srcdir)
subprocess.run(args=[f"cd {srcdir} && make"], shell=True)
sys.path += [srcdir]

