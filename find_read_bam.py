#! /usr/bin/env python
import sys, os, re
import commands

files= commands.getoutput("find ./ -iname \"*tar\"").split()
for f in files:
    l= commands.getoutput("tar -tf "+str(f)).split()
    for item in l:
        if "bai" not in item: 
            print f,item


