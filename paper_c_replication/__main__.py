import argparse
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description='Paper C: downstream replication from precomputed scores.')
    sub=p.add_subparsers(dest='command',required=True)
    demo=sub.add_parser('demo',help='Generate fictional data and run every analysis.')
    demo.add_argument('--work-dir',type=Path,required=True)
    demo.add_argument('--draws',type=int,default=199,help='Demo only. Empirical runs always use 4999.')
    demo.add_argument('--no-figures',action='store_true')
    run=sub.add_parser('run',help='Run on explicitly supplied empirical input files.')
    run.add_argument('--inputs',type=Path,required=True);run.add_argument('--output',type=Path,required=True)
    run.add_argument('--no-figures',action='store_true')
    args=p.parse_args()
    from .analysis import run as execute
    if args.command=='demo':
        from .synthetic import generate
        from .io import fresh_directory
        work=fresh_directory(args.work_dir)
        generate(work/'inputs');execute(work/'inputs',work/'outputs',args.draws,not args.no_figures)
    else:
        execute(args.inputs,args.output,4999,not args.no_figures)

if __name__=='__main__':main()
