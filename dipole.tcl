####################################################################################################
###   Dipole Analyser for VMD :: Probe systems for their dipole orientations and magnitudes.     ###
###   Written in TCL/TK to work with VMD.                                                        ###
###											         ###	
###								                                 ###
###   Funding gratefully provided by Unilever plc and University of Southampton.                 ###
###                                                                                              ###
###   All code is "as is" and there is no implied warranty.                                      ###
###                                                                                              ###
###   Donovan (2012)	                                                                         ###
###                                                                                              ###
####################################################################################################

set z_total 0
set lipids 512
for {set i 1} {$i <= $lipids} {incr i} {
	set z [[atomselect top "name C10 and (resid $i)"] get {z}]
	#puts $z
	set z_total [expr {$z_total+$z}]
}
set midplane [expr {$z_total/$lipids}]
puts "The middle of the bilayer is:  $midplane"

set total_dipole_x 0.0
set total_dipole_y 0.0
set total_dipole_z 0.0
set total_mag_dipole 0.0
set total_mag_cos_theta 0.0
set total_mag_cos_theta_with_lipid 0.0
set count 0
set total_theta 0.0

### Create Hisotgrams

set N 70
for {set i -$N} {$i <= $N} {incr i} {
 set array($i) 0
 set arraylipid($i) 0
 set lipidangle($i) 0 
}
set binwidth [expr {1.0 / $N}]

set n [molinfo top get numframes]
for { set frame 1 } { $frame < $n } { incr frame } {
	puts "Frame: $frame"
	for {set i 1} {$i <= $lipids} {incr i} {
		set z [[atomselect top "name C10 and (resid $i)" frame $frame] get {z}]
		if {$z < $midplane} { 
		  set count [expr {$count+1}]
                  set amide [atomselect top "name H1 O1 C1 O2 and (resid $i)" frame $frame]
		  set dipa [measure dipole $amide -debye]
		  set dipa_x [lindex $dipa 0] 
                  set dipa_y [lindex $dipa 1]  
                  set dipa_z [lindex $dipa 2]
		  set total_dipole_x [expr {$total_dipole_x+$dipa_x}]
		  set total_dipope_y [expr {$total_dipole_y+$dipa_y}]
		  set total_dipole_z [expr {$total_dipole_z+$dipa_z}]
		  set mag_dipa [expr {sqrt($dipa_x*$dipa_x+$dipa_y*$dipa_y+$dipa_z*$dipa_z)}]
		  set total_mag_dipole [expr {$total_mag_dipole + $mag_dipa}]
		  set cos_theta_dipa [expr {($dipa_x * 0.0 + $dipa_y * 0.0 + $dipa_z * 1.0)/$mag_dipa}]
		  set total_mag_cos_theta [expr {$total_mag_cos_theta+$cos_theta_dipa}]

                  ## for the lipid orientation purely. 
                  set h1 [lindex [[atomselect top "name C1 and (resid $i)" frame $frame] get {x y z}] 0]
		  set h2 [lindex [[atomselect top "name C20 and (resid $i)" frame $frame] get {x y z}] 0]
                  set h1x [lindex $h1 0]
		  set h1y [lindex $h1 1]
		  set h1z [lindex $h1 2]
		  set h2x [lindex $h2 0]
		  set h2y [lindex $h2 1]
		  set h2z [lindex $h2 2]
		  set dispvect [vecdist $h1 $h2]
		  set dispnormx [expr {($h1x-$h2x)/$dispvect}] 
                  set dispnormy [expr {($h1y-$h2y)/$dispvect}]
                  set dispnormz [expr {($h1z-$h2z)/$dispvect}]


		  set cos_theta_lipid_angle [expr {($dispnormx * 0.0 + $dispnormy * 0.0 + $dispnormz * 1.0)}] 
                  ## with normal towards middle (check ineq.)
		  set bin [expr {int($cos_theta_lipid_angle/double($binwidth))}]
                  set lipidangle($bin) [ expr { int($lipidangle($bin)) + 1}]
		  set theta [expr {acos($cos_theta_lipid_angle)}]



                  ## for the bilayer normal.
		  set bin [expr {int($cos_theta_dipa/double($binwidth))}]
                  set array($bin) [ expr { int($array($bin)) + 1}]
		  set theta [expr {acos($cos_theta_dipa)}]
                  ## for dot product with lipid director which points towards the middle of the membrane.
                  
		 
		  set cos_theta_dipa_lipid [expr {($dipa_x * $dispnormx + $dipa_y * $dispnormy + $dipa_z * $dispnormz)/$mag_dipa}]
		  set total_mag_cos_theta_with_lipid [expr {$total_mag_cos_theta_with_lipid+$cos_theta_dipa_lipid}]
		  set bin [expr {int($cos_theta_dipa_lipid/double($binwidth))}]
                  set arraylipid($bin) [ expr { int($arraylipid($bin)) + 1}]
		  set theta [expr {acos($cos_theta_dipa_lipid)}]



		  set total_theta [expr {$total_theta+$theta}]
                } else {  }	
	}
}

set av_dipole_x [expr {$total_dipole_x/double($count)}]
set av_dipole_y [expr {$total_dipole_y/double($count)}]
set av_dipole_z [expr {$total_dipole_z/double($count)}]
set mag_dipa_av [expr {sqrt($av_dipole_x*$av_dipole_x+$av_dipole_y*$av_dipole_y+$av_dipole_z*$av_dipole_z)}]
#set cos_theta_dipa_av [expr {($av_dipole_x * 0.0 + $av_dipole_y * 0.0 + $av_dipole_z * 1.0)/double($mag_dipa_av)}]
set mag_dip_av_mag [expr {$total_mag_dipole/$count}]
#set theta_av [expr {acos($cos_theta_dipa_av)}]
set cos_theta_av [expr {$total_mag_cos_theta/$count}]
set theta_av [expr {($total_theta/$count)*(180/3.14159)}]

puts "Dipole |mag| is $mag_dipa_av Debye from super averaging all vectors"
puts "Dipole |mag| is $mag_dip_av_mag Debye from super averaging all magnitudes"
puts "Dipole |cos theta|: $cos_theta_av"
#puts "Direct \theta av: $theta_av"


### Print Histograms

set filename "bilayernormal.txt"
set fileId [open $filename "w"]
for {set i -$N} {$i <= $N} {incr i} {
 set x [ expr { $binwidth*$i} ]
 set freq [expr {$array($i)/double($count)}]
 puts $fileId "$x $freq"
}
close $fileId

set filename "lipiddirector.txt"
set fileId [open $filename "w"]
for {set i -$N} {$i <= $N} {incr i} {
 set x [ expr { $binwidth*$i} ]
 set freq [expr {$arraylipid($i)/double($count)}]
 puts $fileId "$x $freq"
}
close $fileId

set filename "lipidangle.txt"
set fileId [open $filename "w"]
for {set i -$N} {$i <= $N} {incr i} {
 set x [ expr { $binwidth*$i} ]
 set freq [expr {$lipidangle($i)/double($count)}]
 puts $fileId "$x $freq"
}
close $fileId



array unset {}
puts "All done"


	
