import os, sys

def isExist(path, name="File"):
    if os.path.exists(path):
        print("[info] '{}' Exists in '{}'".format(name, path))
    else:
        sys.exit("[error] '{}' Does not exist in path: '{}'".format(name, path))

def isExistOrMake(path, name="File"):
    if os.path.exists(path):
        print("[info] '{}' Exists in '{}'".format(name, path))
    else:
        try:
            print("[info] '{}' Not Exist, creating ...".format(name, path))
            os.mkdir(path)
            print("[info] '{}' Created.".format(name))
        except:
            sys.exit("[error] '{}' Does not exist in path: '{}'".format(name, path))