#!/usr/local/bin/perl

$numargs = @ARGV;
$arg_1 = $ARGV[0];
$arg_2 = $ARGV[1];

if ($numargs < 2) {
   print ("\n");
   print ("  usage: hhvel.pl [input] [ouput]\n");
   print ("         create HASH velocity model\n");
   die   ("\n");
  }

   open(INPUT,"$ARGV[0]");
   open(OUTPUT,">$ARGV[1]");

    while(<INPUT>) {
        ($dum,$vel,$depth) = split(/\s+/,$_);
        if($vel > 0) {
         $depth = sprintf("%04.1f",$depth);
         $vel = sprintf("%6.4f",$vel);
         $model{$depth} = $vel;
         print OUTPUT ("$depth  $vel\n");
         print        ("$depth  $vel\n");
        }
     }  

         foreach $d (keys %model) {
           print " HASH VALUES: $d  $model{$d}\n";
          }

#kds2
# 4.50   0.0
# 5.85   1.0
# 6.00   4.0
# 6.50  10.0
# 6.80  30.0
# 7.85  38.0
