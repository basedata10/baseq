import re
def cleanStr(s):
   s = re.sub('[^0-9a-zA-Z_]', '', s)
   s = re.sub('^[^a-zA-Z_]+', '', s)
   return s