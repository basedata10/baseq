.. _cmd:

命令行
==========
我们使用click库将函数封装成命令。可以在命令行中，在任何路径下对于函数进行调用，方便开发和测试。

Using click
------------
Click 的使用模式如下
::
    import click, os, sys
    CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

    @click.group(context_settings=CONTEXT_SETTINGS)
    @click.pass_context
    def cli(ctx):
        pass

    #define the commands ...
    @cli.command(short_help="tools for CSV files, run 'baseq csv' for help")
    @click.option('--subject', '-s', default='A Email', help='Subject')
    @click.option('--message', '-m', default='The message', help='Message')
    @click.option('--attches', '-a', multiple=True, default='', help='Message')
    def email(subject, message, attches):
        from baseq.utils.email_send import Client
        print(subject, message, attches)
        Client().send_mail(subject, message, attches)

调用命令
::
    baseq email -s "Hello World." -m "Thanks for your help" -a infos.txt