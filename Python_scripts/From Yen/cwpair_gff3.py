# cwpair.py
#
# Watson/Crick strand pairing script
#
# By Pindi Albert, 2011
#
# Input: tab-separated list of called peaks on both strands
# Format: chromosome (chr##), strand (+/-), start (index), end (index), value (read count) 
#
# Output: list of matched pairs and list of unmatched orphans
# Files: S (simple), D (detailed), O (orphans), P (frequency preview plot), F (final frequency plot)
# Files: statistics.txt, summarizes each input file
# Format: see headers
#
# Run with no arguments or -h for usage and command line options

from optparse import OptionParser, IndentedHelpFormatter
import csv, bisect, sys, logging, os, traceback, math

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

try:
    from matplotlib import pyplot
    graph_available = True
    
    # Graph settings, change to customize graph
    Y_LABEL = 'Peak-pair counts'
    X_LABEL = 'Peak-pair distance (bp)'
    TICK_WIDTH = 3
    ADJUST = [0.140, 0.9, 0.9, 0.1] # Amount to shift the graph to make labels fit, [left, right, top, bottom]
    pyplot.rc('xtick.major', size=10.00) # Length of tick marks, use TICK_WIDTH for width
    pyplot.rc('ytick.major', size=10.00)
    pyplot.rc('lines', linewidth=4.00)
    pyplot.rc('axes', linewidth=3.00)
    pyplot.rc('font', family='Arial', size=32.0)
        
except ImportError:
    logging.error('Unable to import matplotlib. Graphs will be unavailable.\n')
    graph_available = False


def distance(peak1, peak2):
    return (peak2[1]+peak2[2])/2 - (peak1[1]+peak1[2])/2
    
    
def gff_row(cname, start, end, score, source, type='.', strand='.', phase='.', attrs={}):
    return (cname, source, type, start, end, score, strand, phase, gff_attrs(attrs))
    
def gff_attrs(d):
    if not d:
        return '.'
    return ';'.join('%s=%s' % item for item in d.items())


def parse_chromosomes(reader):
    chromosomes = {}
    reader.next()
    for line in reader:
        if len(line) == 9: # gff3 format
            cname, junk, junk, start, end, value, strand, junk, junk = line
        else: # txt format
            cname, strand, start, end, value = line
        start = int(start)
        end = int(end)
        value = float(value)
        if cname not in chromosomes:
            chromosomes[cname] = []
        peaks = chromosomes[cname]
        peaks.append((strand, start, end, value))
    return chromosomes
   
def perc95(chromosomes):
    ''' Returns the 95th percentile value of the given chromosomes '''
    values = []
    for peaks in chromosomes.values():
        for peak in peaks:
            values.append(peak[3])
    values.sort()
    return values[int(len(values)*0.95)] # Get 95% value
   
    
def filter(chromosomes, threshold=0.05):
    ''' Filters the peaks to those above a threshold. Threshold < 1.0 is interpreted
    as a proportion of the maximum, >=1.0 as an absolute value.'''
    if threshold < 1:
        p95 = perc95(chromosomes)
        logging.info('The 95th percentile value is %.2f' % p95)
        threshold = p95 * threshold # Make the threshold a proportion of the
    logging.info('The filter cutoff is %.2f' % threshold)
    before = sum([len(values) for values in chromosomes.values()])
    for cname, peaks in chromosomes.items():
        chromosomes[cname] = [peak for peak in peaks if peak[3] > threshold]
    after = sum([len(values) for values in chromosomes.values()])
    logging.info('%d of %d peaks (%d%%) survived filtering' % (after, before, after*100/before))
    
    
def split_strands(chromosome):
    watson = [peak for peak in chromosome if peak[0] == '+']
    crick = [peak for peak in chromosome if peak[0] == '-']
    return watson, crick
    
    
class FrequencyDistribution(object):
    def __init__(self, start, end, binsize=10, d=None):
        self.start = start
        self.end = end
        self.dist = d or {}
        self.binsize = binsize
    def get_bin(self, x):
        ''' Returns the bin in which a data point falls '''
        return self.start + (x-self.start) // self.binsize * self.binsize + self.binsize/2.0
    def add(self, x):
        x = self.get_bin(x)
        self.dist[x] = self.dist.get(x, 0) + 1
    def graph_series(self):
        x = []
        y = []
        for i in range(self.start, self.end, self.binsize):
            center = self.get_bin(i)
            x.append(center)
            y.append(self.dist.get(center, 0))
        return x, y
    def mode(self):
        return max(self.dist.items(), key=lambda data: data[1])[0]
    def size(self):
        return sum(self.dist.values())
        
    
def all_pair_distribution(chromosomes, up_distance, down_distance, binsize):
    dist = FrequencyDistribution(-up_distance, down_distance, binsize=binsize)
    for cname, data in chromosomes.items():
        watson, crick = split_strands(data)
        crick.sort(key=lambda data: float(data[1]))
        keys = make_keys(crick)
        for peak in watson:
            for cpeak in get_window(crick, peak, up_distance, down_distance, keys):
                dist.add(distance(peak, cpeak))
        logging.debug('Processed chromosome %s. Distribution size %d' % (cname, dist.size()))
    return dist


def make_keys(crick):
    return [(data[1] + data[2])//2 for data in crick]
    
    
def get_window(crick, peak, up_distance, down_distance, keys=None):
    ''' Returns a window of all crick peaks within a distance of a watson peak.
    crick strand MUST be sorted by distance'''
    
    strand, start, end, value = peak
    midpoint = (start + end) // 2
    lower = midpoint - up_distance
    upper = midpoint + down_distance
    
    keys = keys or make_keys(crick)
    start_index = bisect.bisect_left(keys, lower)
    end_index = bisect.bisect_right(keys, upper)
    
    return [cpeak for cpeak in crick[start_index:end_index]] 
    


# Match commands functions take a window of crick peaks
#  and the peak to match

def match_largest(window, peak):
    if not window:
        return None
    return max(window, key=lambda cpeak: cpeak[3])
    
def match_closest(window, peak):
    if not window:
        return None
    def key(cpeak):
        d = distance(peak, cpeak)
        if d < 0: # Search negative distances last
            d = 10000 - d # And then prefer less negative distances
        return d
    return min(window, key=key)


def match_mode(window, peak, mode):
    if not window:
        return None
    return min(window, key=lambda cpeak: abs(distance(peak, cpeak)-mode))

METHODS = {'mode':match_mode, 'closest':match_closest, 'largest':match_largest}

PLOT_FORMATS = ['png', 'pdf', 'svg']

colors = 'krg'

# Plot commands

def frequency_plot(freqs, fname, labels=[], title=''):
    if not graph_available:
        logging.warning('Cannot output graph "%s"' % title)
        return
    pyplot.clf()
    pyplot.figure(figsize=(10, 10))
    for i, freq in enumerate(freqs):
        x, y = freq.graph_series()
        pyplot.plot(x, y, '%s-' % colors[i])
    if len(freqs) > 1:
        pyplot.legend(labels)
    pyplot.xlim(freq.start, freq.end)
    pyplot.ylim(ymin=0)
    pyplot.ylabel(Y_LABEL)
    pyplot.xlabel(X_LABEL)
    
    pyplot.subplots_adjust(left=ADJUST[0], right=ADJUST[1], top=ADJUST[2], bottom=ADJUST[3])
    
    ax = pyplot.gca() # get the current axes

    for l in ax.get_xticklines() + ax.get_yticklines():
        l.set_markeredgewidth(TICK_WIDTH)
    
    pyplot.savefig(fname)


def process_file(path, match_func, threshold, options):
    if match_func == 'all':
        funcs = METHODS.keys()
    else:
        funcs = [match_func]
    statistics = [perform_process(path, func, threshold, options) for func in funcs]
    if match_func == 'all':
        frequency_plot([s['dist'] for s in statistics], statistics[0]['graph_path'], labels=METHODS.keys())
    return statistics

def perform_process(path, match_func, threshold, options):
    
    up_distance = options.up_distance
    down_distance = options.down_distance
    binsize = options.binsize
    plot_format = options.plot_format
    
    statistics = {} # Keep track of statistics for the output file
    
    logging.info('Processing file "%s" with match method "%s"' % (path, match_func))

    input = csv.reader(open(path,'rt'), delimiter='\t')
    
    directory, fname = os.path.split(path)
    statistics['fname'] = fname
    statistics['dir'] = directory
    
    if threshold >= 1:
        filter_string = 'fa%d' % threshold
    else:
        filter_string = 'f%d' % (threshold * 100)
        
    attrs = filter_string + 'u%dd%db%d' % (up_distance, down_distance, binsize)
    
    directory = os.path.join(directory, 'cwpair_output_%s_%s' % (match_func, attrs))
    
    if not os.path.exists(directory):
        os.mkdir(directory)
    
    if fname.startswith('INPUT'):
        fname = fname[5:].strip('_') # Strip "INPUT_" from the file if present
    fname = ''.join(fname.split('.')[:-1]) # Strip extension (will be re-added as appropriate)
    
    def make_path(output_type, extension='txt'): # Returns the full path for a certain output
        return os.path.join(directory, '%s_%s.%s' % (output_type, fname, extension))
        
    def td_writer(output_type, extension='txt'): # Returns a tab-delimited writer for a certain output
        return csv.writer(open(make_path(output_type, extension),'wt'), delimiter='\t')
    

    try:
        chromosomes = parse_chromosomes(input)
    except Exception:
        logging.error('Unable to parse file "%s".\n%s' % (path, traceback.format_exc()))
        return

    summary_output = td_writer('S', extension=options.format)
    detailed_output = td_writer('D')
    orphan_output = td_writer('O')
    preview_plot_path = make_path('P', plot_format)
    final_plot_path = make_path('F', plot_format)
    statistics['stats_path'] = os.path.join(directory, 'statistics.txt')
    statistics['graph_path'] = make_path('C', plot_format)
    
    if options.format == 'txt':
        summary_output.writerow(('chrom', 'midpoint', 'midpoint+1','c-w reads sum','c-w distance (bp)'))
    detailed_output.writerow(('chrom', 'start', 'end', 'value', 'strand') * 2 + ('midpoint', 'c-w reads sum', 'c-w distance (bp)'))
    orphan_output.writerow(('chrom', 'strand', 'start', 'end', 'value'))
    
    statistics['perc95'] = perc95(chromosomes)
    if threshold:
        # Apply filter
        filter(chromosomes, threshold)

    
    if match_func == 'mode':
        logging.info('Finding frequency distribution')
        freq =  all_pair_distribution(chromosomes, up_distance, down_distance, binsize)
        mode = freq.mode()
        statistics['preview_mode'] = mode
        logging.info('The mode is %d' % mode)
        logging.info('Outputting frequency preview plot')
        frequency_plot([freq], preview_plot_path, title='Preview frequency plot')
    else:
        statistics['preview_mode'] = 'NA'
            
    logging.info('Pairing peaks')
            
    dist = FrequencyDistribution(-up_distance, down_distance, binsize=binsize)
    orphans = 0
    for cname, chromosome in chromosomes.items():
        logging.debug('Processing chromosome %s' % cname)
            
        # Each peak is (strand, start, end, value)
        watson, crick = split_strands(chromosome)
        logging.debug('%d watson, %s crick peaks' % (len(watson), len(crick)))
    
        watson.sort(key=lambda data: -float(data[3])) # Sort by value of each peak
        crick.sort(key=lambda data: float(data[1])) # Sort by position to facilitate binary search
        
        keys = make_keys(crick)
        
        for peak in watson:
            
            window = get_window(crick, peak, up_distance, down_distance, keys)
            
            if match_func == 'mode':
                match = match_mode(window, peak, mode)
            else:
                match = METHODS[match_func](window, peak)
            
            if match:
                # Write output
                # (chr, start+, end+, value+, +, start-, end-, value-, -, dist, coord)
                midpoint = (match[1] + match[2] + peak[1] + peak[2]) // 4
                d = distance(peak, match)
                dist.add(d)
                if options.format == 'txt':
                    summary_output.writerow((cname, midpoint, midpoint+1, peak[3]+match[3], d))
                else: #gff
                    summary_output.writerow(gff_row(cname=cname, source='cwpair', start=midpoint, end=midpoint+1,
                                                    score=peak[3]+match[3], attrs={'cw_distance':d}))
                detailed_output.writerow((cname,
                                peak[1], peak[2], peak[3], '+', cname, 
                                match[1], match[2], match[3], '-',
                                midpoint, peak[3]+match[3], d))
                i = bisect.bisect_left(keys, (match[1]+match[2])/2)
                del crick[i]
                del keys[i]
            else:
                orphan_output.writerow((cname, peak[0], peak[1], peak[2], peak[3]))
                orphans += 1
               
        for cpeak in crick: # Remaining crick peaks are orphans
            orphan_output.writerow((cname, cpeak[0], cpeak[1], cpeak[2], cpeak[3]))
        orphans += len(crick) 
                
    statistics['paired'] = dist.size() * 2
    statistics['orphans'] = orphans
    logging.info('Paired %d peaks. %d orphans.' % (dist.size() * 2, orphans))
    
    statistics['final_mode'] = dist.mode()
    logging.info('The final mode is %d', dist.mode())
                
    logging.info('Outputting final frequency distribution')
    frequency_plot([dist], final_plot_path, title='Frequency distribution')
    statistics['dist'] = dist
                
    logging.info('Done processing file "%s"' % path)
    
    return statistics
                

def write_statistics(statistics):
    ''' Writes a list of statistics to the file(s) specified by them'''
    logging.info('Writing statistics')
    by_file = {}
    for stats in statistics: # Collect all the stats together by destination file
        if not stats: # Skip "None" statistics from failed files
            continue
        path = stats['stats_path']
        if path not in by_file:
            by_file[path] = []
        by_file[path].append(stats)
                
    for file, statistics in by_file.items():
        f = csv.writer(open(file, 'wt'), delimiter='\t')
        keys = ['fname', 'final_mode', 'preview_mode', 'perc95', 'paired', 'orphans']
        f.writerow(keys)
        for stats in statistics:
            f.writerow([stats[key] for key in keys])
    


usage = '''
input_paths may be:
- a file or list of files to run on
- a directory or list of directories to run on all files in them
- "." to run in the current directory

example usages:
python cwpair.py -m closest /path/to/a/file.txt path/to/another/file.txt
python cwpair.py -u 50 -d 100 -f 5 /path/to/a/data/directory/
python cwpair.py -f 5 -m largest .
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description
 

def run():   
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-m', action='store', type='string', dest='method', default='mode',
                      help='Method of finding match. Valid options are mode (default), closest, largest, and all (run with each method).')
    parser.add_option('-u', action='store', type='int', dest='up_distance', default=50,
                      help='Distance upstream of a Watson peak to allow a Crick pair. Default 50.')
    parser.add_option('-d', action='store', type='int', dest='down_distance', default=100,
                      help='Distance downstream of a Watson peak to allow a Crick pair. Default 100.')
    parser.add_option('-b', action='store', type='int', dest='binsize', default=1,
                      help='Width of bins for frequency plots and mode calculation. Default 1 (no bins).')
    parser.add_option('-f', action='store', type='int', dest='relative_threshold', default=0,
                      help='Percentage of the 95 percentile value to filter below. Default 0 (no filtering).')
    parser.add_option('-F', action='store', type='float', dest='absolute_threshold', default=0,
                      help='Absolute value to filter below. Overrides -f. Default 0 (no filtering).')
    parser.add_option('-p', action='store', type='string', dest='plot_format', default='pdf',
                      help='Format to output graph in. Valid options are png, svg, and pdf (default).')
    parser.add_option('-o', action='store', type='string', dest='format', default='gff',
                      help = 'Format to output simple output file in. Valid options are gff (default) and txt.')
    parser.add_option('-v', action='store_true', dest='verbose', help='Verbose mode: displays debug messages')
    parser.add_option('-q', action='store_true', dest='quiet', help='Quiet mode: suppresses all non-error messages')
    (options, args) = parser.parse_args()
    
    if options.verbose:
        logging.getLogger().setLevel(logging.DEBUG) # Show all info/debug messages
    if options.quiet:
        logging.getLogger().setLevel(logging.ERROR) # Silence all non-error messages
        
    if options.method in METHODS or options.method == 'all':
        match_func = options.method
    else:
        parser.error('%s is not a valid method. Use -h option for a list of valid methods.' % options.method)
        
    if options.plot_format in PLOT_FORMATS:
        plot_format = options.plot_format
    else:
         parser.error('%s is not a valid plot format. Use -h option for a list of valid formats.' % options.plot_format)
         
    if options.absolute_threshold:
        threshold = options.absolute_threshold
    elif options.relative_threshold:
        threshold = options.relative_threshold / 100.0
    else:
        threshold = 0
       
        
    if not args:
        parser.print_help()
        sys.exit(1)
        
    statistics = []
        
    for path in args:
        if not os.path.exists(path):
            parser.error('Path %s does not exist.' % path)
        if os.path.isdir(path):
            for fname in os.listdir(path):
                fpath = os.path.join(path, fname)
                if os.path.isfile(fpath) and not fname.startswith('.'): 
                    stats = process_file(fpath, match_func, threshold, options)
                    statistics.extend(stats)
        else:
            stats = process_file(path, match_func, threshold, options)
            statistics.extend(stats)
            
    write_statistics(statistics)
    
if __name__ == '__main__':
    run()
    #import cProfile
    #cProfile.run('run()', 'profile.bin')
