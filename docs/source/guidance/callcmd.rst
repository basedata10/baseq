软件调用
========
程序大量涉及到调用软件。

我们使用subprocess运行系统命令, 封装成为了 run_cmd函数。调用方式：
::
    from baseq.mgt import run_cmd
    run_cmd("List Folder...", "ls -l")

源代码
::
    from subprocess import call
    import sys

    def run_cmd(name, cmd):
        print("[run] Command: {}".format(cmd))
        try:
            call(cmd, shell=True)
            print("[info] {} complete.".format(name, cmd))
        except:
            sys.exit("[error] {} exit with error.".format(name))
