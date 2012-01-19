set amide [atomselect "name C25 O26 C27 N23 H24 and (resid 1)"]
set n [molinfo top get numframes]
for { set i 0 } { $i < $n } { incr i } {
 measure dipole $amide frame $i
}

set selatoms [atomselect top "resid 2"]

set selatoms [[atomselect top "resid 2 and (resname CER)"] get index]
set selatoms [[atomselect top "name C25 and (resid 5)"] get {x y z}]
set selatoms [[atomselect top "name C25 and (resid 5)"] get {z}]


#### Proper program start

set z_total 0
set lipids 128
for {set i 1} {$i <= $lipids} {incr i} {
	set z [[atomselect top "name C25 and (resid $i)"] get {z}]
	#puts $z
	set z_total [expr {$z_total+$z}]
}
set midplane [expr {$z_total/$lipids}]
puts "The middle of the bilayer is:  $midplane"



set total_amide_x 0.0
set total_amide_y 0.0
set total_amide_z 0.0
set count 0

### Create Hisotgrams

set N 40
for {set i -$N} {$i <= $N} {incr i} {
 set array($i) 0
}
set binwidth [expr {1.0 / $N}]

set n [molinfo top get numframes]
for { set frame 0 } { $frame < $n } { incr frame } {
	puts "Frame: $frame"
	for {set i 1} {$i <= $lipids} {incr i} {
		set z [[atomselect top "name C25 and (resid $i)" frame $frame] get {z}]
		if {$z > $midplane} { 
		  set count [expr {$count+1}]
                  set amide [atomselect top "name C25 C27 N23 H24 O26 and (resid $i)" frame $frame]
                  set glycerol [atomselect top "name C19 C16 C15 C14 C20 O21 H22 and (resid $i)" frame $frame]
		  set dipa [measure dipole $amide -debye]
		  set dipa_x [lindex $dipa 0] 
                  set dipa_y [lindex $dipa 1]  
                  set dipa_z [lindex $dipa 2]
		  set total_amide_x [expr {$total_amide_x+$dipa_x}]
		  set total_amide_y [expr {$total_amide_y+$dipa_y}]
		  set total_amide_z [expr {$total_amide_z+$dipa_z}]


		  set mag_dipa [expr {sqrt($dipa_x*$dipa_x+$dipa_y*$dipa_y+$dipa_z*$dipa_z)}]
		  set cos_theta_dipa [expr {($dipa_x * 0.0 + $dipa_y * 0.0 + $dipa_z * 1.0)/$mag_dipa}]
		  set bin [expr {int($cos_theta_dipa/$binwidth)}]
                  set array($bin) [ expr { int($array($bin)) + 1}]
		  
		  #puts "Amide dipole |magnitude|: $mag_dipa for lipid $i"
		  #puts "Amide dipole cos theta: $cos_theta_dipa for lipid $i"

                } else { #puts "above" }	
	}
}



set av_amide_x [expr {$total_amide_x/double($count)}]
set av_amide_y [expr {$total_amide_y/double($count)}]
set av_amide_z [expr {$total_amide_z/double($count)}]
set mag_dipa_av [expr {sqrt($av_amide_x*$av_amide_x+$av_amide_y*$av_amide_y+$av_amide_z*$av_amide_z)}]
set cos_theta_dipa_av [expr {($av_amide_x * 0.0 + $av_amide_y * 0.0 + $av_amide_z * 1.0)/double($mag_dipa_av)}]
set theta_av [expr {arccos $cos_theta_dipa_av}]

puts "Amide dipole cos theta: $cos_theta_dipa_av at an angle of $theta_av with the vector [1,0,0]"

### Print Hisotgrams



set filename "amide.txt"
set fileId [open $filename "w"]
for {set i -$N} {$i <= $N} {incr i} {
 set x [ expr { $binwidth*$i} ]
 set freq [expr {$array($i)/double($count)}]
 puts $fileId "$x $freq"

 puts "x: $x results: $freq"
}
close $fileId

array unset {}


	
