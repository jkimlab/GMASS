#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use Cwd 'abs_path';
use File::Basename;
use FindBin qw($Bin);
use Getopt::Long;
use Switch;

my $fasta1;
my $fasta2;
my $core = 1;
my $parameter_fle;
my $output_dir = "./";
my $distance = "near";
my @resolutions = ();
my $script_dir = dirname($0);
my $ref_fa; 
my $tar_fa;

GetOptions(
	"f1:s" => \$fasta1,
	"f2:s" => \$fasta2,
	"resolution:s" => \@resolutions,
	"strict:s" => \$distance,
	"param:s" => \$parameter_fle,
	"core:i" => \$core,
	"outdir:s" => \$output_dir,
	"help" => sub{HELP()}
);

my $input_from_cmd = (defined $fasta1 && defined $fasta2 && (scalar @resolutions > 0))? 1 : 0;
my $input_from_paramF = (defined $parameter_fle)? 1 : 0;

if(! $input_from_paramF && ! $input_from_cmd){
	print STDERR "[Error] Check your paramters\n";
	HELP();
}elsif($input_from_cmd && $input_from_paramF){
	print STDERR "[Error] Parameters were received from multiple sources\n";
	HELP();
}

if(defined $parameter_fle){
	open(RF, $parameter_fle);
	while(<RF>){
		chomp;
		if($_ =~ /^\$/){
			my $ln = substr($_, 1);
			my @cols = split(/=/, $ln);
		
			switch($cols[0]){
				case "ASSEMBLY1" {$fasta1 = $cols[1];}
				case "ASSEMBLY2" {$fasta2 = $cols[1];}
				case "RESOLUTION" {@resolutions = split(/,/, $cols[1]);}
				case "STRICT" {$distance = lc($cols[1]);}
				else {
					print STDERR "[Error] Strange paramter: $ln\n";
					exit;
				}
			}
		}
	}
	close RF;
}else{
	@resolutions = split(/,/, join(",", @resolutions));
	$distance = lc($distance);
}

@resolutions = sort {$a <=> $b} @resolutions;
$fasta1 = abs_path($fasta1);
$fasta2 = abs_path($fasta2);
$script_dir = abs_path($script_dir);
$output_dir = abs_path($output_dir);

CHECK_INPUT();

`mkdir -p $output_dir/data`;
# Sef reference and target
my @lengths = ();
my $sequence = "";
my $seq_ID = "";
my $seq_len = 0;
my $total_len = 0;
my $n50_fa1 = 0;
my $n50_fa2 = 0;
my $fa1_newName = "";
my $fa2_newName = "";
my @suffixlist = (".fa", ".fasta");

foreach my $fa ($fasta1, $fasta2){
	my $prefix = basename($fa, @suffixlist);	
	$prefix =~ s/\./_/g;
	`ln -s $fa $output_dir/data/$prefix.fa`;
	`$script_dir/src/kent/faSize -detailed -tab $fa > $output_dir/data/$prefix.size`;

	if($fa eq $fasta1){
		$fa1_newName = $prefix;
	}else{
		$fa2_newName = $prefix;
	}
	
	@lengths = ();
	$sequence = "";
	$seq_len = 0;
	$total_len = 0;

	open(RF, $fa);
	while(<RF>){
		chomp;
	
		if($_ =~ /^>/){
			if($sequence ne ""){
				$seq_len = length($sequence);
				push(@lengths, $seq_len);
				$sequence = "";
				$total_len += $seq_len;
			}
			$seq_ID = substr($_, 1);
		}else{
			$sequence .= $_;
		}
	}
	close RF;
	$seq_len = length($sequence);
	push(@lengths, $seq_len);

	$total_len += $seq_len;
	$total_len /= 2;

	foreach my $len(@lengths){
		$total_len -= $len;
		if($total_len <= 0){
			if($fa eq $fasta1){
				$n50_fa1 = $len;
			}else{
				$n50_fa2 = $len;
			}
			last;
		}
	}
}


if($n50_fa1 >= $n50_fa2){
	$ref_fa = "$output_dir/data/$fa1_newName.fa";
	$tar_fa = "$output_dir/data/$fa2_newName.fa";
}else{
	$ref_fa = "$output_dir/data/$fa2_newName.fa";
	$tar_fa = "$output_dir/data/$fa1_newName.fa";
}

print STDERR "\nSet ". basename ($ref_fa)." as reference, and ". basename($tar_fa) . " as target\n";
print STDERR " -N50 of ". basename ($fasta1) .": $n50_fa1\n";
print STDERR " -N50 of ". basename ($fasta2) .": $n50_fa2\n";

# Detect CSBs
print STDERR "\nDetect CSBs\n";

my $ref_name = basename($ref_fa, @suffixlist);
my $tar_name = basename($tar_fa, @suffixlist);

my $resols = join(",", @resolutions);

my $aln_cmd = "$script_dir/whole_genome_alignment.pl --path $script_dir/path.conf -p $core -res $resolutions[0] -r $ref_fa -t $tar_fa -d $distance -o $output_dir/chainNet 2> $output_dir/aln.log";
my $bld_CSB_cmd = "$script_dir/build_synteny.pl --path $script_dir/path.conf -c $output_dir/chainNet -m $resols -r $ref_name -t $tar_name -o $output_dir/CSB";

print STDERR "$bld_CSB_cmd\n";
if(-d "$output_dir/chainNet"){
	`rm -rf $output_dir/chainNet`;
}

`$aln_cmd`;
`$bld_CSB_cmd`;

# Calculate GMASS
my %stats=();  
my %score = ();
my %scfs = ();
my $result_dir;
my $csb_len = 0;
my $csb_count = 0;
my $scf_len = 0;
my $scf_cnt = 0;
my $used_scf_len = 0;
my $used_scf_cnt = 0;

my $bedtools_path = "";
open(RF, "$script_dir/path.conf");
while(<RF>){
	chomp;

	if($_ =~ /^bedtools=/){
		my @cols = split(/=/);
		$bedtools_path = $cols[1];
		last;
	}
}
close RF;

print STDERR "\nCalculate GMASS\n";
foreach my $res(@resolutions){
	%scfs = ();
	$csb_len = 0;
	$csb_count = 0;
	$scf_len = 0;
	$scf_cnt = 0;
	$used_scf_len = 0;
	$used_scf_cnt = 0;
	
	open(RF, "$output_dir/data/$ref_name.size");
	while(<RF>){
		chomp;
		my @cols = split(/\s+/);
		
		++$scf_cnt;
		my $cur_scf_len = $cols[1];
		$scf_len += $cur_scf_len;
		if($cur_scf_len > $res){
			$used_scf_len += $cur_scf_len;
			++$used_scf_cnt;
		}	
	}

	$result_dir = "$output_dir/CSB/$res/";
	if(-e "$result_dir/Conserved.Segments"){
		chdir $result_dir;
		my $calc_genomecov_ref = "$bedtools_path genomecov -i $ref_name.bed -g $output_dir/data/$ref_name.size > $ref_name.genomecov"; 
		my $calc_genomecov_tar = "$bedtools_path genomecov -i $tar_name.bed -g $output_dir/data/$tar_name.size > $tar_name.genomecov"; 
	
		`$script_dir/syn2bed.pl --path $script_dir/path.conf --sort --merge -s Conserved.Segments`;
		system $calc_genomecov_ref; 
		system $calc_genomecov_tar; 

		$csb_count = `grep ">" Conserved.Segments | wc -l`;
		chomp $csb_count;
	
		my $ln_csb_len = `grep -P "^genome\t1" $ref_name.genomecov`;
		my @cols = split(/\s+/, $ln_csb_len);
		$csb_len = $cols[2];
	}else{
		$csb_count = 0;
		$csb_len = 0;
	}
	
	$stats{"AS_count(ref)"}{$res} = $scf_cnt;
	$stats{"AS_length(ref)"}{$res} = $scf_len;
	$stats{"usedAS_count(ref)"}{$res} = $used_scf_cnt;
	$stats{"usedAS_length(ref)"}{$res} = $used_scf_len;
	$stats{"CSB_length(ref)"}{$res} = $csb_len;
	$stats{"CSB_count(ref)"}{$res} = $csb_count;
	
	%scfs = ();
	$csb_len = 0;
	$scf_len = 0;
	$scf_cnt = 0;
	$used_scf_len = 0;
	$used_scf_cnt = 0;
	
	open(RF, "$output_dir/data/$tar_name.size");
	while(<RF>){
		chomp;
		my @cols = split(/\s+/);
		
		++$scf_cnt;
		my $cur_scf_len = $cols[1];
		$scf_len += $cur_scf_len;
		if($cur_scf_len > $res){
			$used_scf_len += $cur_scf_len;
			++$used_scf_cnt;
		}	
	}
	
	if(-e "$result_dir/Conserved.Segments"){
		my $ln_csb_len = `grep -P "^genome\t1" $tar_name.genomecov`;
		my @cols = split(/\s+/, $ln_csb_len);
		$csb_len = $cols[2];
	}else{
		$csb_len = 0;
	}	

	$stats{"AS_count(tar)"}{$res} = $scf_cnt;
	$stats{"AS_length(tar)"}{$res} = $scf_len;
	$stats{"usedAS_count(tar)"}{$res} = $used_scf_cnt;
	$stats{"usedAS_length(tar)"}{$res} = $used_scf_len;
	$stats{"CSB_length(tar)"}{$res} = $csb_len;
	$stats{"CSB_count(tar)"}{$res} = $csb_count;
}

open (WF, ">$output_dir/assembly_CSB.stats.txt");
my $header = join("\t", @resolutions);
print WF "#stats\t$header\n";
foreach my $st(sort {$a cmp $b} keys %stats){
	print WF "$st";
	foreach my $res(@resolutions){
		if(exists $stats{$st}{$res}){
			print WF "\t$stats{$st}{$res}";
		}else{
			print WF "\tNA";	
		}
	}
	print WF "\n";
}
close WF;

foreach my $res(@resolutions){
	if($stats{"CSB_length(ref)"}{$res} == 0){
		$score{$res}{"Ci score"} = "NA";
		$score{$res}{"Li score"} = "NA";
		$score{$res}{"Si score"} = "NA";
		
		next;
	}

	my $li_score = ($stats{"CSB_length(ref)"}{$res} + $stats{"CSB_length(tar)"}{$res})/($stats{"usedAS_length(ref)"}{$res} + $stats{"usedAS_length(tar)"}{$res});
	
	my $ci_score;
	my $csb_totCnt = $stats{"CSB_count(ref)"}{$res} + $stats{"CSB_count(tar)"}{$res};
	my $asm_totCnt = $stats{"usedAS_count(ref)"}{$res} + $stats{"usedAS_count(tar)"}{$res};
	if($csb_totCnt <= $asm_totCnt){
		$ci_score = 1-((abs($csb_totCnt - $asm_totCnt))/$asm_totCnt);
	}else{
		$ci_score = 1-(abs($csb_totCnt - $asm_totCnt)/((($stats{"usedAS_length(ref)"}{$res}/$res)+($stats{"usedAS_length(tar)"}{$res}/$res))-$asm_totCnt))
	}

	$score{$res}{"Ci score"} = $ci_score;
	$score{$res}{"Li score"} = $li_score;
	$score{$res}{"Si score"} = $ci_score * $li_score;
}

open (WF, ">$output_dir/scores.txt");
print WF "#Resolution\tLi score\tCi score\tSi score\n";
my $si_sum = 0;
my $si_cnt = 0;
foreach my $res(@resolutions){
	print WF "$res\t$score{$res}{'Li score'}\t$score{$res}{'Ci score'}\t$score{$res}{'Si score'}\n";
	if($score{$res}{'Li score'} eq "NA"){
		print STDERR " -The CSBs was not constructed in $res. The GMASS score will be calculated except CSB information in $res.\n";
		next;
	}
	$si_sum += $score{$res}{'Si score'};
	++$si_cnt;
}
print WF "\n--\n";
print WF "GMASS\t". $si_sum/$si_cnt ."\n";
close WF;


##
#subroutines
sub CHECK_INPUT{
	my $pass = 1;
	my %distances = (
			"self" => 1,
			"near" => 1,
			"medium" => 1,
			"far" => 1
	);

	if(! -e $fasta1){
		print STDERR "[Error] Cannot find $fasta1 file\n";
		$pass = 0;
	}
	if(! -e $fasta2){
		print STDERR "[Error] Cannot find $fasta2 file\n";
		$pass = 0;
	}
	if($fasta1 eq $fasta2){
		print STDERR "[Error] The name of fasta file have to be different\n";
		$pass = 0;
	}
	if($#resolutions <= 0){
		print STDERR "[Error] Insufficent resolutions. 2 Resolutions at least are requried\n";
		$pass = 0;
	}
	foreach my $resol(@resolutions){
		if($resol !~ /^\d+?$/){
			print STDERR "[Error] Strange value $resol for resolution. Resolution should be integer value\n";
			$pass = 0;
		}
	}

	if(! exists $distances{$distance}){
		print STDERR "[Error] Strange value $distance for alignment strictness\n";
		$pass = 0;
	}

	if($core == 0){
		print STDERR "[Error] Check the core count\n";
		$pass = 0;
	}

	if($pass){
		my $resols = join(",", @resolutions);
		print STDERR "Input\n";
		print STDERR " -Assembly1: $fasta1\n";
		print STDERR " -Assembly2: $fasta2\n";
		print STDERR " -Resolution: $resols\n";
		print STDERR " -Alignment strictness: $distance\n";
		print STDERR " -cpu: $core\n";
	}else{
		print "\n";
		HELP();

		exit;
	}
}

sub HELP{
	my $src = basename($0);
    print STDERR "Usage: \$ $src [options] -f1 <fasta1> -f2 <fasta2> -r <resolutions> -s <dist>\n";
	print STDERR "\t\t\tor\n";
    print STDERR "       \$ $src [options] -p <param>\n";
    print STDERR "\n -Inputs:\n";
    print STDERR "\t-f1/-f2\t\tUncompressed sequence files in fasta format\n";
    print STDERR "\t-r|--resolution\t\tComma-separated resolution list\n";
    print STDERR "\t-s|--strict\t\tAlignment strictness [self|near|medium|far] (Default: near)\n";
    print STDERR "\t-p|--param\t\tPath of paramter file\n\n";
    
	print STDERR " -Options:\n";
    print STDERR "\t-c|--core\t\tCore number  (Default: 1)\n";
    print STDERR "\t-o|--outdir\t\tPath of output directory  (Default: Current directory)\n";
    print STDERR "\t-h|--help\t\tPrint help message\n";
    exit;
}
