# a simple script to distribute names from a file
import subprocess
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--start', dest="start", action="store",
                    default = [None], type = int, nargs=1,
                    help = "the line of the file (zero-indexed) to start from.")
parser.add_argument('--stop', dest="stop", action="store", default = [None],
                    type = int, nargs=1,
                    help = ("the line of the file (zero-indexed) to stop at. "
                            "(--stop - 1) is the last line processed."))
parser.add_argument('-f','--fname', dest='fname', action="store",
                    required = True, nargs=1,
                    help = ("the file containing a list of names (distributed "
                            "one on each line) that are passed to the script "
                            "as the last line"))

def parse_args():

    if len(sys.argv) <3:
        raise ValueError("need at least 3 argument")

    option_args = [False, False, False]
    
    num_args = 1
    for i, elem in enumerate(sys.argv):
        if i >6:
            num_args=7
            break
        if i % 2 ==1:
            # needs to be a flag
            if elem in ['-f','--fname']:
                if option_args[0]:
                    raise ValueError("can't specify -f/--fname more than once")
                option_args[0] = True
            elif elem == '--start':
                if option_args[1]:
                    raise ValueError("can't specify --start more than once")
                option_args[1] = True
            elif elem == '--stop':
                if option_args[2]:
                    raise ValueError("can't specify --stop more than once")
                option_args[2] = True
            else:
                num_args = i
                break
    else:
        num_args = i+1

    cmd_args = parser.parse_args(sys.argv[1:num_args])
    return cmd_args, sys.argv[num_args:]

def load_list(cmd_args):
    fname = cmd_args.fname[0]
    start = cmd_args.start[0]
    stop = cmd_args.stop[0]

    f = open(fname,'r')
    entries = f.read().splitlines()
    f.close()

    if start is None:
        start = 0
    if stop is None:
        stop = len(entries)

    assert start<stop
    assert 0<=start<len(entries)
    assert stop<=len(entries)
    return entries[start:stop]

def main():
    cmd_args, to_run = parse_args()
    args = load_list(cmd_args)

    if len(to_run)==0:
        for arg in args:
            print arg

    else:
        to_run.append(None)
        for arg in args:
            print "Executing with {:s}".format(arg)
            to_run[-1] = str(arg)
            command = ' '.join(to_run)
            #print command
            subprocess.call(command,shell=True,stdout=sys.stdout,
                            stderr = sys.stderr)
    
if __name__ == '__main__':
    main()
