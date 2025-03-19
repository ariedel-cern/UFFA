import ROOT
import FileUtils as FU

# class that handles the saving of the histograms in the correct file structure
class FileSaver():
    DEBUG = False
    def __init__(self, name = None, option = None):
        self._file = None
        self._wdir  = None
        self._tree  = []
        if name:
            self.touch(name, option)

    # make file and if "new" option is passed, create new file and rename if file already exists
    def touch(self, name, option = None):
        if name[-5:] != '.root':
            name += '.root'
        name = FU.path_expand(name)
        if not option or option == 'new':
            name = FU.file_exists(name)
            self._file = ROOT.TFile(name, 'new')
        else:
            self._file = ROOT.TFile(name, option)
        if not self._file:
            print("FileSaver: file \"" + name + "\" could not be created!")
            return None
        else:
            if FileSaver.DEBUG:
                option = "updated" if option == "updated" else "created"
                print("FileSaver: file \"" + name + "\" " + option + "!")
        self._wdir = self._file
        self._tree.append(self._wdir)

    # write all objects from the list 'objects'
    def write(self, objects):
        for obj in objects:
            obj.Write()

    # make directory in dir_root with name dir_name and save all histos from the list 'objects'
    def writeInDir(self, dir_root, dir_name, objects):
        if not objects:
            return
        dir_new = dir_root.mkdir(dir_name)
        dir_new.cd()
        try:
            self.write(objects)
        except:
            if FileSaver.DEBUG:
                print("\nFileSaver: could not save \"" + objects + "\" in \"" + dir_name + "\"!\n")
            pass
        dir_root.cd()
        del dir_new

    # takes folder input in the form [dir name, [dir content]]
    # where the content can include another folder input
    def write_full_dir(self, folder):
        current = self._wdir
        name = folder[0]
        cont = folder[1]

        if len(name) > 5 and name[-4:] == 'root':
            pass
        else:
            self._wdir = self._wdir.mkdir(name)
            self._wdir.cd()

        for entry in cont:
            if type(entry) == list:
                self.write_full_dir(entry)
            else:
                entry.Write()

        self._wdir = current
        self._wdir.cd()

    # function to change directory
    # cd(0) goes back to the root of the file
    # cd() goes back to the previous dir
    # cd(int) goes back 'int' dir's back
    # cd(str) goes to the 'str' dir
    def cd(self, dir_name = None):
        if not dir_name and dir_name != 0:
            if len(self._tree) == 1:
                return
            else:
                self._tree.pop()                        # remove entry for working directory list
                self._wdir = self._tree[-1]
                return
        else:
            # 0 is the root directory of the file
            if dir_name == 0:
                self._wdir = self._ifile
                self._tree = [self._ifile]              # reset working directory list
                return
            elif type(dir_name) == int:
                for n in range(dir_name):
                    self.cd()
                    if len(self._tree) == 1:
                        return
            else:
                self._set_dir(dir_name)
                return

    # set directory
    def _set_dir(self, dir_name):
        dir_name = FU.path_fix(dir_name)
        name_list = dir_name.rsplit('/')
        for name in name_list:
            dir_new = self._find_obj(name, self._wdir)              # find in current directory
            if not dir_new:
                dir_new = self._find_obj(name, self._ifile)         # find in root directory
                if dir_new:
                    self._tree = [self._ifile]                      # if found in root dir, reset working directory path
            if not dir_new:
                if FileSaver.DEBUG:
                    print("\nDirectory \"" + name + "\" not found!\n")
                    self.ls()
                return None
            self._tree.append(dir_new)                              # append dir to dir path
            self._wdir = dir_new
            if dir_new.Class() == ROOT.TList.Class():
                continue
            dir_new.cd()

    # find object in dir_obj
    def _find_obj(self, obj_name, dir_obj):
        obj = None
        try:                # try to find in TDirectory
            obj = dir_obj.Get(obj_name)
        except:
            pass
        if not obj:         # try to find in tlist
            try:
                obj = dir_obj.findobject(obj_name)
            except:
                if FileSaver.DEBUG:
                    print("\nFileSaver: object \"" + obj_name + "\" not found in \"" + dir_obj + "\"!\n")
                pass
        return obj

    # toggle debug output or set with 'option'
    def setDebug(option = None):
        if type(option) == bool:
            FileSaver.DEBUG = option
            return
        FileSaver.DEBUG = not FileSaver.DEBUG

    # return the file
    def getFile(self):
        return self._file

    # list directory content
    def ls(self):
        self._wdir.ls()

    # list file content
    def fls(self):
        self._file.ls()

    # return working directory
    def pwd(self):
        pwd = ""
        if len(self._tree) == 1:
            return self._file.GetName()
        for obj in self._tree[1:]:
            pwd += obj.GetName() + "/"
        if FileReader.DEBUG:
            print("pwd: \"" + pwd + "\"")
        return pwd

    # return current directory
    def get_dir(self):
        if FileReader.DEBUG:
            print("dir: \"" + self._wdir.GetName() + "\"")
        return self._wdir.GetName()

    def Close(self):
        self._file.Close()
