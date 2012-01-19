

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
set lipids 128
for {set i 1} {$i <= $lipids} {incr i} {
	set z [[atomselect top "name C25 and (resid $i)"] get {z}]
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
set count 0
set total_theta 0.0

### Create Hisotgrams

set N 70
for {set i -$N} {$i <= $N} {incr i} {
 set array($i) 0
}
set binwidth [expr {1.0 / $N}]

set n [molinfo top get numframes]
for { set frame 1 } { $frame < $n } { incr frame } {
	puts "Frame: $frame"
	for {set i 1} {$i <= $lipids} {incr i} {
		set z [[atomselect top "name C25 and (resid $i)" frame $frame] get {z}]
		if {$z > $midplane} { 
		  set count [expr {$count+1}]
                  set amide [atomselect top "name C25 N23 H24 O26 and (resid $i)" frame $frame]
                  set glycerol [atomselect top "name C19 C16 C15 C14 C20 O21 H22 and (resid $i)" frame $frame]
		  set dipa [measure dipole $amide -debye]
		  #puts "Dipole: $dipa"
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
		  set bin [expr {int($cos_theta_dipa/double($binwidth))}]
                  set array($bin) [ expr { int($array($bin)) + 1}]
		  set theta [expr {acos($cos_theta_dipa)}]
		  set total_theta [expr {$total_theta+$theta}]
                } else { #puts "above" }	
	}
}

set av_dipole_x [expr {$total_dipole_x/double($count)}]
set av_dipole_y [expr {$total_dipole_y/double($count)}]
set av_dipole_z [expr {$total_dipole_z/double($count)}]
set mag_dipa_av [expr {sqrt($av_dipole_x*$av_dipole_x+$av_dipole_y*$av_dipole_y+$av_dipole_z*$av_dipole_z)}]
set cos_theta_dipa_av [expr {($av_dipole_x * 0.0 + $av_dipole_y * 0.0 + $av_dipole_z * 1.0)/double($mag_dipa_av)}]
set mag_dip_av_mag [expr {$total_mag_dipole/$count}]
set theta_av [expr {acos($cos_theta_dipa_av)}]
set cos_theta_av [expr {$total_mag_cos_theta/$count}]
set theta_av [expr {($total_theta/$count)*(180/3.14159)}]

puts "Dipole |mag| is $mag_dipa_av Debye from super averaging all vectors"
puts "Dipole |mag| is $mag_dip_av_mag Debye from super averaging all magnitudes"
puts "Dipole |cos theta|: $cos_theta_av"
puts "Direct \theta av: $theta_av"


### Print Hisotgrams

set filename "amide.txt"
set fileId [open $filename "w"]
for {set i -$N} {$i <= $N} {incr i} {
 set x [ expr { $binwidth*$i} ]
 set freq [expr {$array($i)/double($count)}]
 puts $fileId "$x $freq"
}
close $fileId

array unset {}
puts "All done"


	
