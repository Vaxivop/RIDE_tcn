import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--infiles',required=True)
args = parser.parse_args()
print(args.infiles)
a = args.infiles.split(',')
print(type(a))
print(a)
