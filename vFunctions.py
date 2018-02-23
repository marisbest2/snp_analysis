import os
import shutil
import git

class Update_Directory:

    def __init__(self, dependents_dir):
        home = os.path.expanduser("~")

        if os.path.isdir("/bioinfo11/TStuber/Results"):
            #server
            storage_path = "/bioinfo11/TStuber/Results"
            computer_path = "/home/shared"
        elif os.path.isdir("/Volumes/root/TStuber/Results"):
            #mac
            storage_path = "/Volumes/root/TStuber/Results"
            computer_path = "/Users/Shared"

        #check if bioinfo is attached
        if os.path.isdir(storage_path): 
            self.upload_to = storage_path
            self.remote = storage_path + dependents_dir
            self.path = storage_path + dependents_dir
            #make copy of local dependents
            if os.path.isdir(computer_path):
                dep_path = computer_path
                dir_split = dependents_dir.split('/')[1:]
                #build path to file if needed
                for i in dir_split:
                    dep_path += '/' + i
                    if not os.path.exists(dep_path):
                        os.makedirs(dep_path)
                local = computer_path + dependents_dir
                if os.path.isdir(local):
                    try:
                        #remove files
                        shutil.rmtree(local)
                        #copy files from remote to local
                        shutil.copytree(self.remote, local)
                    except:
                        pass
        #if no storage path (bioinfo) check if dependents have been copied local to share
        elif os.path.isdir(computer_path + dependents_dir):
            self.upload_to = "not_found"
            self.remote = "not_found"
            self.path = computer_path + dependents_dir
        #if not in shared folder then check home directory previously downloaded from github
        elif os.path.isdir(home + "/dependencies" + dependents_dir):
            self.upload_to = "not_found"
            self.remote = "not_found"
            self.path = home + "/dependencies" + dependents_dir # sets dependencies directory to home directory
        else:
            os.makedirs(home + "/dependencies")
            print("\n\nDOWNLOADING DEPENDENCIES FROM GITHUB... ***\n\n")
            git.Repo.clone_from("https://github.com/stuber/dependencies.git", home + "/dependencies")
            self.upload_to = "not_found"
            self.remote = "not_found"
            self.path = home + "/dependencies" + dependents_dir