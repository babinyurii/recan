from distutils.core import setup

setup(
  name = 'recan',
  packages = ['recan'],   
  version = '0.1',      
  license='MIT',        
  description = 'recan: recombination analysis tool',   
  author = 'Yuriy Babin',                  
  author_email = 'babin.yurii@gmail.com',      
  url = 'https://github.com/babinyurii/recan', 
  download_url = 'https://github.com/babinyurii/recan/archive/v_0.1.tar.gz',
  keywords = ['DNA recombination', 'bioinformatics', 'genetic distance'],   
  install_requires=[            
          'pandas',
          'plotly',
          'biopython',
          'matplotlib'
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      
    'Intended Audience :: Developers',      
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   
    'Programming Language :: Python :: 3',      
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
)
