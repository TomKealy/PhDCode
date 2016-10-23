Documentation of Slasso(Structured Lasso) Matlab package
by Rodolphe Jenatton 
For any questions or report any bugs, please contact me at rodolphe [DOT] jenatton [AT] inria [DOT] fr

DESCRIPTION

This package is a Matlab program that implements methods for structured variable selection, based on sparsity-inducing norms. 
For more information, please read the following paper:

(2009) R. Jenatton, J.-Y. Audibert and F. Bach.  Structured variable selection with sparsity-inducing norms. Technical report, arXiv:0904.3523

INSTALLATION

You will need first to install the free SOCP solver SDPT3 (available at http://www.math.nus.edu.sg/~mattohkc/sdpt3.html).

Then, specify in SETUP_TOOLBOX.m  the path where the SOCP solver SDPT3 has been installed.

At the Matlab prompt, type "SETUP_TOOLBOX" and all mex files will be compiled. 
Detailed demo scripts are:

   + Demo_for_groups.m
   + Demo_for_slasso.m
   + Demo_for_islasso.m
   + Demo_for_activeset.m

This version has been tested with R2007b,R2008a Matlab versions, on 64-bit Linux machines.

Note that if you want to reproduce the Lasso experiments of the paper, you will need the free implementation of Lasso developed by J.Mairal 
(available at http://www.di.ens.fr/willow/SPAMS/).


LICENSE

The Slasso Matlab package is provided free for non-commercial use under the terms of the GNU General Public License. 
The current version is the 1.0 and has been released on September, 16th 2009. 

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
