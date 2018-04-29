import base64, json

class pack_web_datas():
    def __init__(self):
        self.data = {}

    def add_image(self, name, path):
        print("[info] Add image {}:{}".format(name, path))
        try:
            with open(path, "rb") as image_file:
                encoded_string = base64.b64encode(image_file.read())
                self.data[name] = encoded_string.decode('ascii')
        except:
            print("[error] Failed to get the image file")

    def write_file(self, outpath):
        print("[info] The result file is write to '{}'".format(outpath))
        try:
           with open(outpath, "w") as file:
               json.dump(self.data, file)
        except:
            print("[info] Failed to open file to write, '{}'".format(outpath))