#!/usr/local/bin/perl

$numargs = @ARGV;
$arg_1 = $ARGV[0];
$arg_2 = $ARGV[1]; if ($numargs < 3) {
   print ("\n");
   print ("  usage: model_interp.pl [input] [Limits] [Interp]\n");
   print ("         create HASH velocity model\n");
   die   ("\n");
  }

   open(INPUT,"$ARGV[0]");
   open(OUTPUT,">$ARGV[1]");
   open(INTERP,">$ARGV[2]");

    while(<INPUT>) {
        ($dum,$vel,$depth) = split(/\s+/,$_);

        if($depth > -0.1 && $depth < 0.1) {
          $dz = sprintf("%5.2f",$depth);
          $vz = sprintf("%6.4f",$vel);
          print INTERP "$dz   $vz\n";
        }
         
        if($vel > 0) {
 #       $depth = sprintf("%04.1f",$depth);
 #       $vel = sprintf("%5.2f",$vel);
         push(@v,$vel);
         push(@z,$depth);
        }
     }

      $z_ref = \@z;
      $v_ref = \@v;
      for ($i=1; $i <= $#$z_ref; $i++) {
       $vdiff = $v[$i] - $v[$i-1];
       $zdiff = $z_ref->[$i] - $z_ref->[$i-1];
       print "$v[i]  $z[$i]  ZDIF: $zdiff VDIF: $vdiff\n";
       $inc = $vdiff/$zdiff;
#      $zdiff = sprintf("%5.2f",$zdiff);
#      $vdiff = sprintf("%5.2f",$vdiff);
#      $inc = sprintf("%4.2f",$inc);

       print OUTPUT "ZD: $zdiff  IC: $inc VD: $vdiff  Z: $z[$i-1] $z[$i]  V: $v[$i-1] $v[$i]\n";

       $ddepth   = $z[$i-1];
       $velnew  =  $v[$i-1];
#      print INTERP "$velnew   $ddepth  $inc\n";

        while ($ddepth < $z[$i])  {
          $velnew   = $velnew + $inc;
          $ddepth  = $ddepth + 1;
          $dz = sprintf("%5.2f",$ddepth);
          $vz = sprintf("%6.4f",$velnew);
          print INTERP "$dz   $vz\n";
         }

     }

#kds2
# 4.50   0.0
# 5.85   1.0
# 6.00   4.0
# 6.50  10.0
# 6.80  30.0
# 7.85  38.0
#ZD:  1.00  IC: 1.35 VD:  1.35  Z: 00.0 01.0  V: 4.50 5.85
#ZD:  3.00  IC: 0.05 VD:  0.15  Z: 01.0 04.0  V: 5.85 6.00
#ZD:  6.00  IC: 0.08 VD:  0.50  Z: 04.0 10.0  V: 6.00 6.50
#ZD: 20.00  IC: 0.01 VD:  0.30  Z: 10.0 30.0  V: 6.50 6.80
#ZD:  8.00  IC: 0.13 VD:  1.05  Z: 30.0 38.0  V: 6.80 7.85
#ZD: 22.00  IC: 0.02 VD:  0.35  Z: 38.0 60.0  V: 7.85 8.20
