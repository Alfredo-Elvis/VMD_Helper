import math,os
#last update 2025-02-28

xtc_file = "Your trajectory file[xtc format]"
pdb_file = "You structure file [pdb format]"
model_ID = -1 #start from 0
Show_Residues_Near_MOL = True

def analysis_pdb(file):
    ATOMS=[]
    chains=[]
    residues=["ASH","CLY","LYS","ARG","HIS","HIE","HIP","HID","ASP","GLU","GLH","SER","THR","ASN","GLN","CYS","CYM","SEC","GLY","PRO","ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP","NME","ACE"]
    solution=["SOL","HOH","CA","NA","CL","ZN"]
    chain=["A","B","C","D","E","F","G","H","I","J","K","L","M","N"]
    with open(file,"r") as f:
        lines=f.readlines()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                ATOMS.append(line)
    flag=True
    i=0
    atom_index = -1
    lines=[]
    chain_name=chain[i]
    for atom in ATOMS:
        atom_index += 1
        residue_name=atom[17:20].strip()
        atom_name=atom[12:16].strip()
        if residue_name in residues:
            if atom_name=="CA":
                if flag:
                    flag=False
                    chain_name=chain[i]
                    [x1,y1,z1]=[float(atom[30:38]),float(atom[38:46]),float(atom[46:54])]
                else:
                    chain_name=chain[i]
                    [x,y,z]=[float(atom[30:38]),float(atom[38:46]),float(atom[46:54])]
                    dis=math.sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)
                    [x1,y1,z1]=[x,y,z]
                    if dis>5:
                        i+=1
                        chain_name=chain[i]
                        residue_number=atom[22:26].strip()
                        temp_index = -15
                        for temp in ATOMS[atom_index+temp_index:atom_index]:
                            residue_number_temp=temp[22:26].strip()
                            if residue_number_temp==residue_number:
                                lines[temp_index]=temp[:21]+chain_name+temp[22:]
                            temp_index += 1
        elif residue_name in solution:
            chain_name="S"
        else:
            chain_name="X"
        if chain_name not in chains:
            chains.append(chain_name)
        atom=atom[:21]+chain_name+atom[22:]
        lines.append(atom)
    chain_name_default="A"
    f=open("clean.pdb","w")
    for line in lines:
        chain_name=line[21]
        if chain_name != chain_name_default:
            chain_name_default=chain_name
            f.write("TER\n")
        f.write(line)
    f.write("TER\nEND\n")
    f.close()
    return chains,chain

def analysis_hbond(chain,time="all"):
    current_working_directory = os.getcwd().replace("\\","/")
    if "X" in chain:
        sel1 = ""
        sel2 = "chain X"
        for c in chain:
            if "X" not in c:
                sel1= sel1 + f"chain {c} or "
        command = f"hbonds -sel1 [atomselect {model_ID} \"{sel1[:-4]}\"] -sel2 [atomselect {model_ID} \"{sel2}\"] -writefile yes -upsel yes -frames {time} -dist 3.5 -ang 40 -plot no -log {current_working_directory}/hbonds.log -outfile {current_working_directory}/hbonds.dat -polar yes -DA both -type unique -detailout {current_working_directory}/hbonds_tailout.dat"
        return command
    else:
        return None

def chain_paint(chains,chain,model_ID):
    global Show_Residues_Near_MOL
    i=0
    lines=[]
    color = 0
    for c in chains:
        color_id = 31 - color
        color = color + 1
        if c in chain:
            #Protein
            line=f'''mol addrep {model_ID}
mol modselect {i} {model_ID} chain {c}
mol modcolor {i} {model_ID} ColorID {color_id}
mol modmaterial {i} {model_ID} AOChalky
mol modstyle {i} {model_ID} NewCartoon 0.22 50.00 3.87 1
mol smoothrep {model_ID} {i} 3

'''
            if "X" in chains and Show_Residues_Near_MOL:
                line=line + f'''mol addrep {model_ID}
mol modselect {i+1} {model_ID} (same residue as (protein exwithin 5 of chain X)) and chain {c} and noh
mol modmaterial {i+1} {model_ID} AOChalky
mol modstyle {i+1} {model_ID} CPK 1.40 0.70 50.00 50.00
mol smoothrep {model_ID} {i+1} 3

mol addrep {model_ID}
mol modselect {i+2} {model_ID} ((same residue as (protein exwithin 5 of chain X)) and chain {c}) and name \\"C.*\\"
mol modcolor {i+2} {model_ID} ColorID {color_id}
mol modmaterial {i+2} {model_ID} AOChalky
mol modstyle {i+2} {model_ID} CPK 1.40 0.70 50.00 50.00
mol smoothrep {model_ID} {i+2} 3

'''
                i = i + 2
        elif c == "S":
            #Solution
            line=f'''mol addrep {model_ID}
mol modselect {i} {model_ID} chain {c}
mol modcolor {i} {model_ID} Name
mol modmaterial {i} {model_ID} AOChalky
mol modstyle {i} {model_ID} VDW 0.60 50.00
mol smoothrep {model_ID} {i} 3

'''
        elif c == "X":
            #MOL
            line=f'''mol addrep {model_ID}
mol modselect {i} {model_ID} chain {c} and noh
mol modcolor {i} {model_ID} Name
mol modmaterial {i} {model_ID} AOChalky
mol modstyle {i} {model_ID} CPK 1.40 0.70 50.00 50.00
mol smoothrep {model_ID} {i} 3

mol addrep {model_ID}
mol modselect {i+1} {model_ID} chain {c} and name \\"C.*\\"
mol modcolor {i+1} {model_ID} ColorID {32}
mol modmaterial {i+1} {model_ID} AOChalky
mol modstyle {i+1} {model_ID} CPK 1.40 0.70 50.00 50.00
mol smoothrep {model_ID} {i+1} 3

'''
            i = i + 1
        lines.append(line)
        i = i + 1
    return lines

def painting(xtc_file,pdb_file,model_ID=1):
    pathway = os.getcwd().replace("\\","\\\\")
    (chains,chain)=analysis_pdb(pdb_file)
    lines=chain_paint(chains,chain,model_ID)
    default=f'''color Display Background white
axes location Off

color change rgb 26 0.929412 0.415686 0.352941
color change rgb 27 0.576471 0.862745 0.898039
color change rgb 28 0.980392 0.752941 0.368627
color change rgb 29 0.709804 0.568627 0.937255
color change rgb 30 1.0 0.85098 0.447059
color change rgb 31 0.866667 0.890196 0.843137
color change rgb 32 0.619608 0.898039 0.482353

cd {pathway}

display resetview
mol new {{clean.pdb}} type {{pdb}} first 0 last -1 step 1 waitfor -1
mol addfile {{{xtc_file}}} type {{xtc}} first 0 last -1 step 1 waitfor -1 {model_ID}
mol rename {model_ID} {{{pdb_file[:-4]}}}
mol delrep 0 {model_ID}

'''
    f=open(f"input.tcl","w")
    f.write(default)
    for line in lines:
        f.write(line)
    f.write("display resetview\n\n")
    hbond = analysis_hbond(chains)
    if hbond != None:
        f.write(hbond)
    f.close()
    return True

if __name__ == "__main__":
    if model_ID == -1:
        print("Please set the correct model_ID, start from 0.")
        exit(0)
    if xtc_file == "Your trajectory file[xtc format]" or xtc_file.split(".")[-1] != "xtc":
        print("Please specify your trajectory file[xtc format].")
        exit(0)
    if pdb_file == "You structure file [pdb format]" or pdb_file.split(".")[-1] != "pdb":
        print("Please specify your structure file[pdb format].")
        exit(0)
    try:
        flag = painting(xtc_file,pdb_file,model_ID)
    except Exception as e:
        flag = False
        print(f"Error: {e}")
    if flag:
        print("Please input 'source input.tcl' in the tk Console")
    else:
        print("Please adjust input parameters or input molecules manually.")
