#!/usr/bin/perl
#usage: perl script.pl input.fasta
#purpose: read a fasta file for unique haplotypes/isoforms, cluster exact sequences, include exact partial matches

open(IN, '<', @ARGV[0]);

$fasta = do{local $/; <IN>};
@fasta = split(/\>/, $fasta);

%uniq = undef;
%table = undef;
%counts = undef;

#build up hash of uniq sequences and store seqIDs
foreach $fasrec (@fasta) {
	unless ($fasrec eq undef) {
		@data = split(/\n/, $fasrec, 2);
		@temp = split(/\s+/, @data[0]);
		@data[0] = @temp[0];
		@data[1] =~ s/\-|\n|\t|\s//sgi;
		unless (@data[1] eq undef) {
			if ($uniq{@data[1]} eq undef) {
				$uniq{@data[1]} = @data[0];
				###print @data[1] . " " . $uniq{@data[1]} . "\n";
			}
			else {
				$uniq{@data[1]} .= " " . @data[0];
				###print @data[0] . " " . @data[1] . "\n";
				###print @data[1] . " " . $uniq{@data[1]} . "\n";
			}
		}
		###print $fasrec . "\n";
	}
}

=begin - this step was a bad idea - metaeuk does not indicate truncations, which will get combined by this method... removed this feature (for now)
#check hash of uniq sequences for exact overlaps 
foreach $rec1 (keys %uniq) {
	foreach $rec2 (keys %uniq) {
		unless ($rec1 eq undef || $rec2 eq undef) {
			if ($rec1 =~ /$rec2/ && $rec1 ne $rec2) {
				$uniq{$rec2} .= " " . $uniq{$rec1};
				###print $rec1 . " " . $rec2 . "\n\n";
				delete $uniq{$rec1};
			}
			elsif ($rec2 =~ /$rec1/ && $rec2 ne $rec1) {
				$uniq{$rec1} .= " " . $uniq{$rec2};
				delete $uniq{$rec2};
			}
			print "";
		}
	}
}
=cut

#output1: uniqseq / ids
open(OUT1, '>', @ARGV[0] . ".byuniqseq");
foreach $rec1 (sort {$a cmp $b} keys %uniq) {
	unless ($rec1 eq undef) {
		print OUT1 $rec1 . "\t" . $uniq{$rec1} . "\n";
	}
}

#output2: list of id / uniqseq
open(OUT2, '>', @ARGV[0] . ".byseqid");

foreach $rec1 (sort {$a cmp $b} keys %uniq) { unless ($rec1 eq undef) {
	@temp = split(/\s+/, $uniq{$rec1});
	foreach $rec2 (@temp) { unless ($rec2 eq undef) {
		print OUT2 $rec2 . "\t" . $rec1 . "\n";
		@temp2 = split(/\|/, $rec2);
		$table{@temp2[0]}->{$rec1} += 1;
		$counts{$rec1}++;
	}}
}}

#output3: table of uniqseq x id
open(OUT3, '>', @ARGV[0] . ".byseqid_table");
print OUT3 "ID";
#foreach $rec1 (sort {$a cmp $b} keys %uniq) { unless ($rec1 eq undef) {
foreach $rec1 (sort {$counts{$b} <=> $counts{$a}} keys %counts) { unless ($rec1 eq undef) {
	print OUT3 "\t$rec1";
}}
print OUT3 "\n";
foreach $id (sort {$a cmp $b} keys %table) { unless ($id eq undef) {
	print OUT3 "$id";
	#foreach $rec1 (sort {$a cmp $b} keys %uniq) { unless ($rec1 eq undef) {
	foreach $rec1 (sort {$counts{$b} <=> $counts{$a}} keys %counts) { unless ($rec1 eq undef) {
		if ($table{$id}->{$rec1} eq undef) {
			$table{$id}->{$rec1} = 0;
		}
		print OUT3 "\t" . $table{$id}->{$rec1};
	}}
	print OUT3 "\n";
}}


exit;