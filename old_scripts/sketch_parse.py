#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('sketch_json', type=argparse.FileType('r'))
args = parser.parse_args()

import json 
d = json.load(args.sketch_json)
print(d)
