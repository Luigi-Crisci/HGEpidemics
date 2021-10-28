import os

def backup_files(path):
	for file in os.scandir(path):
		if file.is_dir():
			backup_files(file)
			continue
		dest_path = path.path + "/old_results/"
		if not os.path.exists(dest_path):
			os.mkdir(dest_path)
		os.rename(file, dest_path + file.name)