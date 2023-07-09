import subprocess
import sys
import os

srcdir = __file__.split('/')[:-1] + "/src"
subprocess.run([f"cd {srcdir};", "make"], shell=True)
sys.path.append(srcdir)

