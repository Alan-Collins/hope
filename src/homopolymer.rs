#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]

#[derive(Debug)]
pub struct HomopolymerRecord {
    pub contig: String,
    pub start: u32,
    pub stop: u32,
    pub base: String,
    pub length: u32,
}

impl HomopolymerRecord {
    pub fn print(&self) {
        println!("{:?}", self);
    }
}