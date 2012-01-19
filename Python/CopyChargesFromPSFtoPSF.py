####################################################################################################
###   CopyChargesFromPStoPSF :: Looks up charges and creates a dictionary which enables transfer ### 
###   of charge information between file formats.                                                ###
###											         ###	
###								                                 ###
###   Funding gratefully provided by Unilever plc and University of Southampton.                 ###
###                                                                                              ###
###   All code is "as is" and there is no implied warranty.                                      ###
###                                                                                              ###
###   Donovan (2011)  Flagged up for a rewrite....                                               ###
###                                                                                              ###
####################################################################################################

import cfg
mydict = {}

def Get_Charges():
	FILE_frame = open(cfg.input_file, "r")
	line = FILE_frame.readline().strip()
	
	started_bool = False
	quit = False
	while (quit == False):   
		if (started_bool) and (line == ""):
			quit = True
			break
			#print "breaking..."
		#print line
		if (line != ""):
			index = str(line.split()[0])
			#print index
			if (index) == "1":
				started_bool = True

                	 
		if (started_bool) and (line != ""):
			chem = str(line.split()[4])
			charge = str(line.split()[6])
			#print chem + " " + charge
			mydict[chem] = charge
		line = FILE_frame.readline().strip()	
	FILE_frame.close()
	print mydict

def Replace_Charges():
	new_segment = ""
	FILE_frame = open(cfg.output_file, "r")
	FILE_write = open(cfg.write_file, "w")
	line = FILE_frame.readline().strip()
	started_bool = False
	quit = False
	while (quit == False):   
		if (started_bool) and (line == ""):
			quit = True
			break
			print "breaking..."
		if (line != ""):
			index = str(line.split()[0])
			if (index) == "1":
				started_bool = True

		if (started_bool) and (line != ""):
			chem = str(line.split()[4])
			charge = str(line.split()[6])
			A_curr = int(line.split()[0])
			B_curr = str(line.split()[1])
			C_curr = str(line.split()[2])
			D_curr = str(line.split()[3])
			E_curr = str(line.split()[4])
			F_curr = str(line.split()[5])
			G_curr = float(line.split()[6])
			H_curr = float(line.split()[7])
			I_curr = str(line.split()[8])
			if chem in mydict:
				
				G_curr = float(mydict[chem])
				
			new_line = "   %5d" % A_curr + " " + B_curr + "   " + C_curr.ljust(5," ") + "" + D_curr.ljust(3," ") + "  " + E_curr.ljust(3," ") + "  " + F_curr + "   %+5f" % G_curr + "       %+5f" % H_curr + "         " + I_curr + "\n"
			
		
			FILE_write.write(new_line)
		line = FILE_frame.readline().strip()	
		
	FILE_write.close()
	FILE_frame.close()


Get_Charges()
Replace_Charges()
