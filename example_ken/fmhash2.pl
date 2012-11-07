: # use perl
eval 'exec $ANTELOPE/bin/perl -S $0 "$@"'
if 0;

use lib "$ENV{ANTELOPE}/data/perl" ;
use Datascope;
use Math::Trig;
use List::Util 'max';
use List::Util 'min';
use POSIX qw(log10);

$numargs = @ARGV;

if ($numargs < 4) {
   print ("\n");
        print "Usage: fmhash.pl [Year] [Jday] [Orid] [Output File]\n";
   die   ("\n");
  }

#Example '>fm.pl 2008 114 507701'     M 2.9 near Mogul 2008

   $year = $ARGV[0];
   $jday = $ARGV[1];
   $evtorid = $ARGV[2];
   $outputfile= $ARGV[3];

   open(OUTPUT,">$ARGV[3]");

   $jday = sprintf("%3.3d",$jday);
   $DB = "/data\/$year\/$jday/reno";

      if(-e $DB) {

# Join tables in dbprocess

      @db  = dbopen("$DB", "r");

      @dbday = dbprocess ( @db, "dbopen origin",
                 "dbsubset orid == $evtorid",
                 "dbjoin assoc",
                 "dbjoin arrival",
                 "dbjoin site",
                 "dbjoin affiliation",
                 "dbsubset iphase =~ /.*[PS].*/");

# Count all phases
#     @dbsub = dbsubset( @dbday , "fm =~ /.*[cd].*/");

#      $fms = dbquery( @dbsub , "dbRECORD_COUNT" );   
       $fms = dbquery( @dbday , "dbRECORD_COUNT" );   
      print " $DB $fms $evtorid\n";

      $ifms = 0; 

# Event summary information. 

      while ($ifms < $fms) {
        $dbday[3] = $ifms; 
        ($net,$sta,$chan,$iphase,$fm,$ordtime,$lat,$lon,$ml,$depth) = dbgetv(@dbday,"affiliation.net","arrival.sta","arrival.chan","arrival.iphase",
                "arrival.fm","origin.time","origin.lat","origin.lon","origin.ml","origin.depth");
         $fm = substr($fm, 0, 1);

         if($fm eq 'c') {
           $fm = "U";
          }
         elsif($fm eq 'd') {
           $fm = "D";  
         }

         &format_text;

         if($ifms < 1 && $chan =~ /..Z/) {
           print "$phasecard\n";
           print "$sta $net  $chan I $fm\n";
           print OUTPUT "$phasecard\n";
           print OUTPUT "$sta $net  $chan I $fm\n";
          }
         elsif ($chan =~ /..Z/) {
           print        "$sta $net  $chan I $fm\n";
           print OUTPUT "$sta $net  $chan I $fm\n";
          }

         $ifms++;
       }
        dbclose (@db);
      }  
#    }

#--------------------------------------------------------------

   sub format_text {

         ($latd,$latm) = unpack("A2 A5",$lat);
         ($lond,$lonm) = unpack("x A3 A5",$lon);

#        print "Source Location:  $latd $latm $lond $lonm\n";

         $latm = $latm * 60;
         $lonm = $lonm * 60;
         $latm = sprintf("%5.2f",$latm);
         $lonm = sprintf("%5.2f",$lonm);

#        print "Source Loc Format:  $latd $latm $lond $lonm\n";

         $depth = sprintf("%05.2f",$depth);
         $timeres = sprintf("%7.4f",$timeres);
         $sta = sprintf("%-4.4s",$sta);
         $ordtime  = `epoch +"%Y%m%e%H%M%S.%s" $ordtime`;
         chop($ordtime);
#        $phtime  = sprintf("%15.2f",$phtime);
#        $ordtime = sprintf("%15.2f",$ordtime);
         $phasecard = sprintf("%15.2f%2.2d %5.2f%3.3d %5.2f%05.2f    %2.2d",$ordtime,$latd,$latm,$lond,$lonm,$depth,$fms);
       }

#1994 12111 415.5034 14.55118 37.0618.13    31
