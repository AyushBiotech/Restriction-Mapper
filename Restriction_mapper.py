import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Restriction import AllEnzymes
from Bio.Restriction import RestrictionBatch
import pandas as pd
from datetime import datetime
print("="*67)
print(" [ RESTRICTION MAPPER ] Scanning Workspace....")
print("="*67)

print("[1] Input DNA from file")
print("[2] Type / paste DNA manually")

choice = input("Choose input method: ").strip()

if choice == "1":
    print("[1] Android")
    print("[2] Windows")
    print("[3] Linux")
    print("[4] macOS")

    choice = input("Select your operating system: ").strip()

    if choice == "1":
        base = "/storage/emulated/0"
        print(" Nice to see you my fellow Android user😁")

    elif choice == "2":
        base = os.path.expanduser("~")

    elif choice == "3":
        base = os.path.expanduser("~")

    elif choice == "4":
        base = os.path.expanduser("~")

    else:
        print("Invalid choice, try running it again ")
        exit()

    print(f" Your base directory is : [ {base} ]")

    folder_name = input("Enter your folder path relative to the base directory : ").strip()

    folder_path = os.path.join(base, folder_name)

    while True:
        try:
            folder = os.listdir(f"{folder_path}")
            break

        except FileNotFoundError:
            print(" This Folder Doesn't exist, kindly check your input")
            folder_name = str(input(" Kindly enter your Folder path again : ")).strip()
            folder_path = os.path.join(base, folder_name)

        except PermissionError:
            print(" Sorry!! You don't have permission to access this folder ")
            folder_name = str(input(" Kindly enter your Folder path again : ")).strip()
            folder_path = os.path.join(base, folder_name)
    allowed = [".txt",".fasta",".fa",".gb",".gbk",".genbank"]

    if folder == []:
        print(" This folder is empty, kindly put your input files in the folder")
        exit()

    else:
        allowed_files = []

        for items in folder:
            name, extension = os.path.splitext(items)

            if extension.lower() in allowed:
                allowed_files.append(items)

        if allowed_files == []:
            print(" No supported file in this folder 📂")
            exit()

        for index, allowed_items in enumerate(allowed_files, start=1):
            print(f" [{index}] {allowed_items}")

        while True:
            try:
                choice = int(input(" Enter your choice :"))
                real_input = choice - 1
                file_path = allowed_files[real_input]
                break

            except IndexError:
                print(" Invalid choice, please select a listed file  ")

            except ValueError:
                print(" Invalid choice, Choice must be an integer ")

        path = os.path.join(folder_path, file_path)

        print(f" The Full file path is :[ {path} ] ")
        name, extension = os.path.splitext(path)
        extension = extension.lower()
        if extension in [".fasta", ".fa"]:
        	try :
        		record = SeqIO.read(path, "fasta")
        		dna = record.seq
        	except ValueError :
        		print(" Selected File is Empty or Invalid")
        		exit()
        elif extension in [".gb", ".gbk", ".genbank"]:
        	try:
        		record = SeqIO.read(path, "genbank")
        		dna = record.seq
        	except ValueError :
        		print(" Selected File is Empty or Invalid")
        		exit()
        elif extension == ".txt":
        	with open(path,"r") as file:
        		raw_data = file.read()
        	raw_data = raw_data.strip().upper().replace(" ", "").replace("\n","").replace("\r","")
        	if raw_data == "" :
        		print(" Selected File is Empty or Invalid")
        		exit()
        	dna = Seq(raw_data)

elif choice == "2":
    dna = Seq(input("Paste or type your DNA sequence : ").strip().upper().replace(" ", "").replace("\n", ""))
    if dna == "":
    	print(" No Input is given, cancelling the analysis... ")
    	exit()

else:
    print("Invalid choice, try running it again ")
    exit()
print("="*67)
print(" [Analysing The Input]....")
print("="*67)
bases =["A","T","G","C"]
for base in dna :
	if base not in bases :
		print(" Invalid DNA is entered!!! ")
		exit()
print(" DNA Topology :")
print(" [1] Linear")
print(" [2] Circular")
topology = (input(" Choose Topology of the DNA :")).strip()
if topology == "1":
	shape = True
elif topology == "2":
	shape = False
else :
	print(" Invalid Choice")
	exit()
timestamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
res_batch = []
length = len(dna) 
print(f" \n•The length of the DNA sequence is : {length}")
while True :
		enzyme_name = input(" Enter a Restriction Enzyme OR type All to check for every possible enzyme( OR press Enter to finish) : ").strip()
		if enzyme_name == "":
			break
		if enzyme_name == "All" or enzyme_name == "all" or enzyme_name == "ALL":
			res_batch =(AllEnzymes)
			break
		elif enzyme_name in AllEnzymes:
			res_batch.append(enzyme_name)
		else :
			print(" No Enzyme with this name Exists!! Please Try again ")
if len(res_batch) > 0:
	my_batch = RestrictionBatch(res_batch)
	analysis = my_batch.search(dna,linear = shape)
	rows = []
	for enzyme, position in analysis.items():
		if len(position) == 0 :
			print(f" \n[{enzyme}] :\n Total Sites found = 0\n Non-cutter")
			rows.append([enzyme, len(position), "None", "Non-cutter"])
		if len(position) != 0 and enzyme.is_blunt() == True:
			print(f" \n[{enzyme}] :\n Total Sites found = {len(position)}\n Cuts the sequence at {position}\n Cut Type = Blunt")
			rows.append([enzyme, len(position), position, "Blunt"])
		if len(position) != 0 and enzyme.is_3overhang() == True:
			print(f" \n[{enzyme}] :\n Total Sites found = {len(position)}\n Cuts the sequence at {position}\n Cut Type = 3' Overhang ")
			rows.append([enzyme, len(position), position, "3' Overhang"])
		if len(position) != 0 and enzyme.is_5overhang() == True:
			print(f" \n[{enzyme}] :\n Total Sites found = {len(position)}\n Cuts the sequence at {position}\n Cut Type = 5' Overhang ")
			rows.append([enzyme, len(position), position, "5' Overhang"])	
	print("="*67) 
	my_input = str(input (" Do  you  want  to  filter  the  analysis  on the basis of cutting frequency(Yes or No only) : ")).lower().strip()
	print("="*67)
	if my_input == "yes":
		while True :
			try:
				cut_frequency = int(input (" Enter your desired Cutting Frequency :"))
				break
			except ValueError :
				print("Cutting Frequency must be an  integer!!! Please try again")
		found = False
		for enzyme, position in analysis.items() :
			if len(position) == (cut_frequency) and enzyme.is_blunt() == True:
				found = True
				print(f" \n[{enzyme}] :\n Total Sites found = {len(position)}\n Cuts the sequence in {position}\n Cut Type = Blunt")
			if len(position) == (cut_frequency) and enzyme.is_3overhang() == True:
				found = True
				print(f" \n[{enzyme}] :\n Total Sites found = {len(position)}\n Cuts the sequence in {position}\n Cut Type = 3' Overhang")
			if len(position) ==(cut_frequency) and enzyme.is_5overhang()== True :
				found = True
				print(f" \n[{enzyme}] :\n Total Sites found = {len(position)}\n Cuts the sequence in {position}\n Cut Type = 5' Overhang")
		if found == False:
			print(" \n No enzyme with this cutting frequency exists for this sequence")
	
	elif my_input == "no":
		pass
	else :
		print(" Choice must be Yes Or No ")
	print("="*67)
	input2 = str(input("Do you want a CSV analysis report ?(Yes or No only) : ")).lower().strip()
	print("="*67)
	if input2 == "yes" :
		df = pd.DataFrame(rows, columns=[
    "Enzyme",
    "Total Sites",
    "Cutting Positions",
    "Cut Type"
])
		df.to_csv(f"restriction_analysis({timestamp}).csv", index=False)
		print(" This full  analysis  report is saved in a CSV formated with a timestamp")
	elif input2 == "no" :
		exit()
	else :
		print(" Choice must be Yes Or No ")		
else :
	print(" [No Enzyme name is entered, analysis is cancelled ]")
