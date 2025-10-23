// process SEQTK_ORIENT {
//     fair true
//     tag "$query"
//     label 'process_low'
//     publishDir(
//       path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
//       mode: 'copy',
//       overwrite: true,
//       saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
//     ) 
//     input:
//         tuple val(reference), val(query), path(query_genome), path(alignment_info)

//     output:
//         tuple val(query), path("*oriented.fa"), emit: oriented

//     script:
//         """
//         cat ${alignment_info} | awk '{ if (\$9 == -1) {print \$11}}' > inverted_seqs.txt
//         cat ${alignment_info} | awk '{ if (\$9 == 1) {print \$11}}' > forward_seqs.txt
//         seqtk subseq ${query_genome} inverted_seqs.txt > ${query}_rev.fa
//         seqtk subseq ${query_genome} forward_seqs.txt > ${query}_fwd.fa
//         seqtk seq -r ${query}_rev.fa > ${query}_rev_rc.fa
//         cat ${query}_fwd.fa ${query}_rev_rc.fa > ${query}_oriented.fa
//         """
// }
// process SEQTK_SUBSET {
//     fair true
//     tag "$meta"
//     label 'process_low'
//     publishDir(
//       path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
//       mode: 'copy',
//       overwrite: true,
//       saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
//     ) 

//     input:
//         tuple val(meta), path(genome)

//     output:
//         tuple val(meta), path("*_subset.fa"), emit: subset

//     script:
//     def pattern = params.subset_pattern
//         """
//         if [[ ${genome} == *.gz ]];
//             then
//                 gzip -dc ${genome} > ${meta}_genome.fa
//         else
//             ln -s ${genome} ${meta}_genome.fa
//         fi
//         grep $pattern ${meta}_genome.fa | sed 's/>//' > names.lst
//         seqtk subseq ${meta}_genome.fa names.lst > ${meta}_subset.fa
//         """
// }

process SEQTK_ORIENT {
    fair true
    tag "$query"
    label 'process_low'
    publishDir(
        path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() },
        mode: 'copy',
        overwrite: true,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    )

    input:
    tuple val(reference), val(query), path(query_genome), path(alignment_info)

    output:
    tuple val(query), path("*oriented.fa"), emit: oriented

    script:
    """
    set -euo pipefail

    # columns assumed: strand in \$9 (1 or -1); seq name in \$11
    awk '{ if (\$9 == -1) print \$11 }' ${alignment_info} > inverted_seqs.txt
    awk '{ if (\$9 ==  1) print \$11 }' ${alignment_info} > forward_seqs.txt

    # If a list is empty, create an empty fasta so concatenation still works
    if [[ -s inverted_seqs.txt ]]; then
      seqtk subseq ${query_genome} inverted_seqs.txt > ${query}_rev.fa
      seqtk seq -r ${query}_rev.fa > ${query}_rev_rc.fa
    else
      : > ${query}_rev_rc.fa
    fi

    if [[ -s forward_seqs.txt ]]; then
      seqtk subseq ${query_genome} forward_seqs.txt > ${query}_fwd.fa
    else
      : > ${query}_fwd.fa
    fi

    cat ${query}_fwd.fa ${query}_rev_rc.fa > ${query}_oriented.fa
    """
}
process SEQTK_SUBSET {
  tag "$name"
  label 'process_low'

  publishDir(
    path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() },
    mode: 'copy', overwrite: true,
    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
  )

  input:
    tuple val(name), path(fasta)

  output:
    tuple val(name), path("${name}_subset_*.fa")

  // capture the pattern; default to '.*'
  script:
  def pattern  = (params.subset_pattern?.toString()?.trim() ?: '.*')
  def patt_md5 = java.security.MessageDigest.getInstance('MD5').digest(pattern.bytes).encodeHex().toString()
  def patt_hash = patt_md5.substring(0, 8)   // stable 8-char hash

  """
  set -euo pipefail

  pattern=${groovyShellQuote(pattern)}
  patt_hash=${groovyShellQuote(patt_hash)}

  # Extract IDs matching the regex
  grep '^>' "${fasta}" | sed 's/^>//; s/ .*//' | grep -E "\$pattern" > ids.txt

  # Safety: if ids.txt is empty, fail early with a clear message
  if [ ! -s ids.txt ]; then
    echo "ERROR: No headers matched subset_pattern: \$pattern" >&2
    exit 1
  fi

  # Subset
  seqtk subseq "${fasta}" ids.txt > "${name}_subset_\${patt_hash}.fa"
  """
}

// Helper to safely quote strings into bash
def groovyShellQuote(s){ return "'" + s.replace("'", "'\"'\"'") + "'" }