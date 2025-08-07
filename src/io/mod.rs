/*
This module provides the I/O functions (e.g. VCF/BCF readers, text readers etc).
*/
use anyhow::Result;
use rust_htslib::bcf;
use std::path::Path;
use std::{
    fmt::Display,
    io::{BufRead, BufReader, Read},
};

// Define multople readers for the indexed and unindexed XCF file
pub enum XcfReader {
    Indexed(bcf::IndexedReader),
    Readthrough(bcf::Reader),
}

/// Same as smakcr, but single threaded for now
pub fn read_xcf<P: AsRef<Path> + Display>(path: P, has_index: bool) -> Result<XcfReader> {
    let xcf_reader: XcfReader = if has_index {
        XcfReader::Indexed(
            bcf::IndexedReader::from_path(path).expect("Cannot load indexed BCF/VCF file"),
        )
    } else {
        XcfReader::Readthrough(bcf::Reader::from_path(path).expect("Cannot load BCF/VCF file"))
    };
    Ok(xcf_reader)
}

/// Same as smakcr, but single threaded for now
pub fn read_file<P: AsRef<Path> + Display>(
    path: P,
) -> Result<impl Iterator<Item = Result<String>>> {
    let sniffed_reader: std::result::Result<(Box<dyn Read>, niffler::Format), niffler::Error> =
        niffler::from_path(&path);
    // if the file has fewer than 5 bytes, `niffler` can't sniff the compression format and will
    // return a `FileTooShort` error; this could be due to
    // * an empty file
    // * a file containing only a single FASTA record with the ID consisting only of a single
    //   character and the sequence being empty
    // * a file containing a single sequence with a one-character ID and one-character sequence and
    //   missing newline character at the end
    // we don't want to fail at this stage in these cases and thus handle the `FileTooShort` error
    // separately
    let reader = match sniffed_reader {
        Ok(rdr) => Ok(rdr.0),
        Err(e) => Err(e),
    }?;
    let buffer = BufReader::new(reader);
    let records = buffer
        .lines()
        .map(|item: std::result::Result<String, std::io::Error>| {
            item.map_err(|e| anyhow::anyhow!("Error reading line: {}", e))
        });
    Ok(records)
}
