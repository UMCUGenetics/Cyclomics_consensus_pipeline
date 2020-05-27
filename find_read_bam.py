#! /usr/bin/env python
import subprocess

files = subprocess.getoutput("find ./ -iname \"*tar\"").split()
for f in files:
    l = subprocess.getoutput("tar -tf {f}".format(f=f)).split()
    for item in l:
        if "bai" not in item: 
            print(f, item)


