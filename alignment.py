'''
Created on Mar 18, 2016

@author: covingto
'''


import os, os.path, sys

class AlignmentInterface(object):
    ##
    # returns the source sequence in source sequence space (- included)
    def get_source_sequence(self, start = None, end = None):
        raise NotImplementedError()
    
    ##
    # returns the query sequence in source sequence space (- included)
    def get_query_sequence(self, start = None, end = None):
        raise NotImplementedError
    
    @staticmethod
    def conform_bases(b1, b2):
        assert len(b1) == len(b2), "Base strings are not the same length"
        b1_new = []
        b2_new = []
        for i in range(len(b1)):
            if b1[i] != '-' or b2[i] != '-':
                b1_new.append(b1[i])
                b2_new.append(b2[i])
        return ''.join(b1_new), ''.join(b2_new)
    
class BaseAlignment(AlignmentInterface):
    def __init__(self, scontig, sstart, ssize, sstrand, ssrcSize, sbases, qcontig, qstart, qsize, qstrand, qsrcSize, qbases):
        self.scontig = scontig
        self.sstart = int(sstart)
        self.ssize = int(ssize)
        self.sstrand = sstrand
        self.ssrcSize = int(ssrcSize)
        self.qcontig = qcontig
        self.qstart = int(qstart)
        self.qsize = int(qsize)
        self.qstrand = qstrand
        self.qsrcSize = int(qsrcSize)
        
        sbases, qbases = self.conform_bases(sbases, qbases)
        
        self.sbases = sbases 
        self.qbases = qbases 
    
    def get_source_sequence(self, start=None, end=None):
        if start is None:
            start = self.sstart
        if end is None:
            end = self.sstart + self.ssize
        if start == self.sstart and end == self.sstart + self.ssize:
            return self.sbases
        else:
            sindex, eindex = self.get_offset_indexes(start, end)
            return self.sbases[sindex:eindex]
    
    def get_query_sequence(self, start=None, end=None):
        if start is None:
            start = self.sstart
        if end is None:
            end = self.sstart + self.ssize
        if start == self.sstart and end == self.sstart + self.ssize:
            return self.qbases
        else:
            sindex, eindex = self.get_offset_indexes(start, end)
            return self.qbases[sindex:eindex]
    
    def get_offset_indexes(self, start, end):
        assert start < end, "Start %s is not less than end %s" % (str(start), str(end))
        assert start >= self.sstart, "Start %s is outside of this alignment %s" % (str(start), str(self.sstart))
        assert end <= self.sstart + self.ssize, "End is outside of this alignment"
        # remember that - don't consume bases so our size must account for this 
        eindex = sindex = None
        counter = self.sstart
        for i, v in enumerate(self.sbases):
            if counter == start:
                sindex = i
            if v != '-':
                counter += 1
            if counter > end:
                eindex = i
                break
        return sindex, eindex
    
    def substart(self):
        return self.sstart
    
    def subend(self):
        return self.sstart + self.ssize
    
class CompoundAlignment(AlignmentInterface):
    def __init__(self, alignments):
        self.alignments = alignments
    
    def get_query_sequence(self, start=None, end=None):
        if start is None:
            start = self.alignments[0].substart()
        if end is None:
            end = self.alignments[-1].subend()
        seq = ''
        for a in self.alignments:
            if a.subend() < start: continue 
            if a.substart() > end: break 
            _t_start = max(start, a.substart())
            _t_end = min(end, a.subend())
            seq += a.get_query_sequence(_t_start, _t_end)
        return seq
    
    def get_source_sequence(self, start=None, end=None):
        if start is None:
            start = self.alignments[0].substart()
        if end is None:
            end = self.alignments[-1].subend()
        seq = ''
        for a in self.alignments:
            if a.subend() < start: continue 
            if a.substart() > end: break 
            _t_start = max(start, a.substart())
            _t_end = min(end, a.subend())
            seq += a.get_source_sequence(_t_start, _t_end)
        return seq
    
    def substart(self):
        return self.alignments[0].substart()
    
    def subend(self):
        return self.alignments[-1].subend()

class AlignmentBuilder(object):
    def __init__(self):
        self._alignments = []
    
    def add(self, align):
        self._alignments.append(align)
    
    def build(self):
        return CompoundAlignment(self._alignments)
    
class AlignmentContainer(object):
    def __init__(self, reader, subject = 'hg19', query = 'mm10'):
        self._alignments = []
        self._load_alignments(reader, subject, query)
        
    
    def _load_alignments(self, reader, subject, query):
        next_line = reader.readline()
        s_match = ('scontig', 'sstart', 'ssize', 'sstrand', 'ssrcSize', 'sbases')
        q_match = ('qcontig', 'qstart', 'qsize', 'qstrand', 'qsrcSize', 'qbases')
        my_block = {}
        while next_line != '': # not the same as '\n'
            lsplit = [c.strip() for c in next_line.split()]
            if len(lsplit) < 1 or lsplit[0] == '':
                # new paragraph
                if my_block:
                    if 'qcontig' not in my_block:
                        my_block.update(dict(zip(q_match, ['none', 0, 0, my_block.get('sstrand'), 0, '-' * len(my_block.get('sbases'))])))
                    self.make_alignment(my_block)
                    my_block = {}
            elif lsplit[0] == 's':
                if subject in lsplit[1]:
                    my_block.update(dict(zip(s_match, lsplit[1:])))
                elif query in lsplit[1]:
                    my_block.update(dict(zip(q_match, lsplit[1:])))
            elif lsplit[0] == 'e' and query in lsplit[1]:
                if 'sbases' in my_block:
                    sbases = my_block.get('sbases')
                    eclass = lsplit[-1]
                    if eclass == 'C':
                        lsplit[-1] = '-' * len(sbases)
                    else:
                        lsplit[-1] = '=' * len(sbases)
                    my_block.update(dict(zip(q_match, lsplit[1:])))
            next_line = reader.readline()
    
    def make_alignment(self, block):
        if len(block) == 12:
            self._alignments.append(BaseAlignment(**block))
        else:
            print block
            
    ##
    # return an alignment object representing the local multiseq alignment
    # for this region in the subject
    #
    # @param start is the zero based start location of the alignment
    # @param end is the zero based end location of the alignment
    def get_alignment(self, start, end):
        # first we need the index of the start in the _alignments
        start_index = self._binary_start_search(start)
        if start_index is None:
            return None
        alignment_builder = AlignmentBuilder()
        _index = start_index
        while True:
            _align = self._alignments[_index]
            alignment_builder.add(_align)
            if end > _align.subend():
                _index += 1
                continue
            else:
                break
        return alignment_builder.build()
    
    def _binary_start_search(self, position, start = 0, end = None):
        if end is None:
            end = len(self._alignments)
        if end < start:
            return None
        mid = (start + end) / 2
        _align = self._alignments[mid]
        if _align.substart() > position:
            return self._binary_start_search(position, start, mid - 1)
        elif _align.subend() < position:
            return self._binary_start_search(position, mid + 1, end)
        else:
            return mid
