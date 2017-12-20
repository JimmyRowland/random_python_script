import os

def rename(path):
    for subdir in get_immediate_subdirectories(path):
        for filename in os.listdir(subdir):
            # print(filename[-4:])
            # print(path+"/"+filename)
            # print(os.path.isfile(path+"/"+filename))
            # print(os.path.isdir(path+"/"+filename))
            if ".mp4" != filename[-4:] and os.path.isfile(os.path.join(subdir,filename)):
                os.rename(os.path.join(subdir,filename), subdir+"/"+""+filename+".mp4")
                print(path + "/" + filename)


# rename("/media/toor/60E2-1EFC/video/videos/may 8")
def get_immediate_subdirectories(a_dir):
    return [os.path.join(a_dir, name) for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

rename("/media/toor/storage/ubuntu download/video")