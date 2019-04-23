#!/usr/bin/perl

use strict;
use warnings;
use FindBin;
use Cwd;
use Cwd 'abs_path';
use File::Basename;
use Parallel::ForkManager;
use Getopt::Long qw(:config no_ignore_case);
use FindBin '$Bin';
use Switch;

my $chainNet_dir;
my $resolution;
my $ref_name;
my $tar_name;
my $path_conf;
my $out_dir;

GetOptions (
	"path|t=s" => \$path_conf,
	"res|m=s" => \$resolution,
	"chainNet|c=s" => \$chainNet_dir,
	"ref|r=s" => \$ref_name,
	"tar|t=s" => \$tar_name,
	"outdir|o=s" => \$out_dir,
);

$chainNet_dir = abs_path($chainNet_dir);

my $makeblocks_src = "";
open (PATH,$path_conf);
while (<PATH>){
	chomp;
	next if /^#/;
	next if /""/;
	my ($program, $path) = split (/=/, $_);
	switch ($program){
		case("makeBlocks"){   $makeblocks_src = $path;  }
	}
}
close(PATH);

`mkdir -p $out_dir`;
my @resolutions = split(/,/,$resolution);

foreach my $res (@resolutions){
	`mkdir -p $out_dir/$res`;
	`sed -e 's:<willbechanged>:$makeblocks_src:' $makeblocks_src/DATA/Makefile > $out_dir/$res/Makefile`;
	open(CONFIG,">$out_dir/$res/config.file");
	print CONFIG ">netdir\n$chainNet_dir\n";
	print CONFIG ">chaindir\n$chainNet_dir\n\n";
	print CONFIG ">species\n";
	print CONFIG "$ref_name\t0\t0\n";
	print CONFIG "$tar_name\t1\t0\n";
	print CONFIG "\n>resolution\n$res\n";
	close(CONFIG);
	`make -C $out_dir/$res`;
	`make tidy -C $out_dir/$res`;
}
