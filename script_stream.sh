#!/bin/sh

# List of supported species and their folder names
species_names="homo_sapiens 
			   mus_musculus 
			   drosophila_melanogaster 
			   saccharomyces_cerevisiae 
			   aquila_chrysaetos_chrysaetos
			   anas_zonorhyncha
			   electrophorus_electricus
			   loxodonta_africana
			   callorhinchus_milii
			   dromaius_novaehollandiae
			   bubo_bubo
			   sciurus_vulgaris
			   accipiter_nisus
			   dicentrarchus_labrax
			   gadus_morhua"
species_paths="homo_sapiens 
			   mus_musculus 
			   drosophila_melanogaster 
			   saccharomyces_cerevisiae 
			   aquila_chrysaetos_chrysaetos
			   anas_zonorhyncha
			   electrophorus_electricus
			   loxodonta_africana
			   callorhinchus_milii
			   dromaius_novaehollandiae
			   bubo_bubo
			   sciurus_vulgaris
			   accipiter_nisus
			   dicentrarchus_labrax
			   gadus_morhua"

# common process steps
process_stream(){
  # Remove FASTA header lines (starting with '>')
  grep -v '^>' | \
  # Convert to uppercase
  tr '[:lower:]' '[:upper:]' | \
  # Keep only ACTG characters
  tr -Cd ACTG
}

print_help(){
  echo "Usage: $0 [-s <species_name>] [-l <limit_bytes>]"
  echo ""
  echo "Description:"
  echo "  Streams and cleans a genome sequence from a URL or standard input."
  echo "  Use -s to fetch from a URL, or pipe data via stdin if -s is omitted."
  echo ""
  echo "Input Options:"
  echo "  -s   Species name (fetches from Ensembl URL)"
  echo ""
  echo "Other Options:"
  echo "  -l   Limit in bytes for URL downloads (optional)"
  echo "  -h   Show this help message"
  echo ""
  echo "Available species:"
  for sp in $species_names; do echo "  - $sp"; done
}

# Parse args
species=""
limit=""
while getopts "s:l:h" opt; do
  case $opt in
    s) species=$OPTARG ;;
    l) limit=$OPTARG ;;
    h) print_help; exit 0 ;;
    *) print_help; exit 1 ;;
  esac
done

# case 1: Read from a species URL (-s)
if [ -n "$species" ]; then
  # Verify species is supported
  found="false"
  i=1
  for name in $species_names; do
    if [ "$species" = "$name" ]; then
      found="true"
      j=1
      for path in $species_paths; do
        [ $j -eq $i ] && species_path=$path && break
        j=$((j+1))
      done
      break
    fi
    i=$((i+1))
  done

  [ "$found" = "false" ] && echo "Unsupported species: $species" >&2 && print_help && exit 1

  base_url="https://ftp.ensembl.org/pub/release-110/fasta/$species_path/dna/"

  # Try to find a valid FASTA file
  file_name=$(curl -s "$base_url" | grep -oE '[^"]+\.dna\.primary_assembly\.fa\.gz' | head -n1)

  # Fallback to toplevel
  [ -z "$file_name" ] && file_name=$(curl -s "$base_url" | grep -oE '[^"]+\.dna\.toplevel\.fa\.gz' | head -n1)

  [ -z "$file_name" ] && echo "No suitable FASTA file found for $species" >&2 && exit 1

  full_url="$base_url$file_name"

  # Build and run curl + filters
  curl_opts="-sL"
  [ -n "$limit" ] && curl_opts="$curl_opts --range 0-$limit"

  curl $curl_opts "$full_url" 2>/dev/null \
  | gunzip -c 2>/dev/null \
  | process_stream

# case 2: Read from stdin
else
  # if no -s is given, assume data is being piped in
  process_stream
fi