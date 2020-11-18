import os
import re
import sys

def select_file(dpath, key = r".*"): # keyは正規表現
    dflist = os.listdir(dpath)
    flist = [f for f in dflist if os.path.isfile(os.path.join(dpath, f)) and re.match(key, f)]
    while(True):
        print("q: quit")
        for i, fname in enumerate(flist):
            print(str(i)+": "+fname)
        select = input("Select file number: ")
        if select.lower() == "q":
            raise ValueError("quit")
        try:
            file_idx =  int(select)
            if(0<=file_idx and file_idx < len(flist)):
                return flist[file_idx]
            else:
                continue
        except ValueError:
            continue
        

def main():
    try:
        #select_file("./config")
        fname = select_file("./config", r".*\.json")
        #select_file("./config", r".*\.csv")
        print(fname)
    except ValueError as q:
        print(q.args[0])
        return
    return

if __name__ == "__main__":
    main()