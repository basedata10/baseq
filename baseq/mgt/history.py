from baseq.setting import history_command_file
import readline, os

def history_cmd():
    print("History Commands...")
    print(history_command_file)

def add_record(dir, cmd=""):
    #print(history_command_file)
    cmd = readline.get_history_item(0)
    length = readline.get_current_history_length()
    cwd = os.getcwd()
    print(length, cwd, cmd)