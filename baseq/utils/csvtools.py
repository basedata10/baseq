

doc = """
Tools for Manipulate CSV Files...
#slice columns of csv...
baseq csv slicecol 4,5 ./example.csv 

"""

def csvtools(command):
    if len(command) == 0:
        print(doc)
    else:
        print(command)