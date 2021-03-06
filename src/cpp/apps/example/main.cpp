#include <iostream>
#include <string>

#include <random>

#include "image/Image.h"
#include "image/Image_funcs.h"
#include "image/Image_storage.h"

#include "DiffeoFunctionMatching.h"

using ImageLib::TImage;

void run_and_save_example(TImage<double>* I0, TImage<double>* I1, const std::string& subpath, const std::string& desc)
{
  printf("%s Initializing\n", subpath.c_str());

  DiffeoFunctionMatching::create(nullptr, nullptr, 0, 0, 0, false);

  bool compute_phi = true;
  std::unique_ptr<DiffeoFunctionMatching> dfm;
  std::string msg;
  std::tie(dfm, msg) = DiffeoFunctionMatching::create(I0, I1, 0.1, 0.1, 0.1, compute_phi);

  /*
  dm = difforma_base.DiffeoFunctionMatching(
    source=I0, target=I1,
    alpha=0.001, beta=0.03, sigma=0.05
  )*/

  printf("%s Running\n", subpath.c_str());
  //dm.run(1000, epsilon=0.1)
}


/*
def run_and_save_example(I0, I1, subpath, description):
  """Utility function to run and export results for a test case."""
  print('"%s": Initializing' % subpath)
  dm = difforma_base.DiffeoFunctionMatching(
    source=I0, target=I1,
    alpha=0.001, beta=0.03, sigma=0.05
  )
  print('"%s": Running' % subpath)
  dm.run(1000, epsilon=0.1)
  
  print('"%s": Creating plots' % subpath)
  if not path.exists(subpath):
    os.makedirs(subpath)
  
  # 2x2 overview plot
  plt1 = plt.figure(1, figsize=(11.7,9))
  plt.clf()

  plt.subplot(2,2,1)
  plt.imshow(dm.target, cmap='bone', vmin=dm.I0.min(), vmax=dm.I0.max())
  plt.colorbar()
  plt.title('Target image')

  plt.subplot(2,2,2)
  plt.imshow(dm.source, cmap='bone', vmin=dm.I0.min(), vmax=dm.I0.max())
  plt.colorbar()
  plt.title('Template image')

  plt.subplot(2,2,3)
  plt.imshow(dm.I, cmap='bone', vmin=dm.I0.min(), vmax=dm.I0.max())
  plt.colorbar()
  plt.title('Warped image')

  plt.subplot(2,2,4)
  # Forward warp    
  phix = dm.phix
  phiy = dm.phiy
  # Uncomment for backward warp
  #phix = dm.phiinvx
  #phiy = dm.phiinvy
  plot_warp(phix, phiy, downsample=4)
  plt.axis('equal')
  warplim = [phix.min(), phix.max(), phiy.min(), phiy.max()]
  warplim[0] = min(warplim[0], warplim[2])
  warplim[2] = warplim[0]
  warplim[1] = max(warplim[1], warplim[3])
  warplim[3] = warplim[1]

  plt.axis(warplim)
  plt.gca().invert_yaxis()
  plt.gca().set_aspect('equal')
  plt.title('Warp')
  plt.grid()
  plt1.savefig(path.join(subpath, 'overview.png'), dpi=300, bbox_inches='tight')

  # Energy convergence plot
  plt2 = plt.figure(2, figsize=(8,4.5))
  plt.clf()
  plt.plot(dm.E)
  plt.grid()
  plt.ylabel('Energy')
  plt2.savefig(os.path.join(subpath, 'convergence.png'), dpi=150, bbox_inches='tight')

  # Dedicated warp plot (forward only)
  plt3 = plt.figure(3, figsize=(10,10))
  plt.clf()
  plot_warp(phix, phiy, downsample=4, )
  plt.axis('equal')
  warplim = [phix.min(), phix.max(), phiy.min(), phiy.max()]
  warplim[0] = min(warplim[0], warplim[2])
  warplim[2] = warplim[0]
  warplim[1] = max(warplim[1], warplim[3])
  warplim[3] = warplim[1]

  plt.axis(warplim)
  plt.gca().invert_yaxis()
  plt.gca().set_aspect('equal')
  plt.title('Warp')
  plt.axis('off')
  #plt.grid(color='black')
  plt3.savefig(path.join(subpath, 'warp.png'), dpi=150, bbox_inches='tight')
  
  # Output description
  with open(path.join(subpath, 'description.txt'), 'w') as f:
    f.write(subpath)
    f.write('\n')
    f.write(description)
    
  print('Done at ' + time.asctime() + '\n')
*/


void test1() {
  //std::random rnd;
  //std::mt19937_64 gen(rnd());
  std::mt19937_64 gen(1234);
  std::string description =
  "Translations ought to be exactly achievable even with periodic\n"
  "boundary conditions. This test verifies that presumption.\n\n"
  "It seems that images that are non-smooth on the pixel level causes divergence.\n"
  "Some binary \"smoothing\" method may be employed. Instead of single pixels,\n"
  "small squares or circles could be used.";
  int nPoints = 30;
  int delta = 20;
  auto I0 = ImageLib::create_filled_image<double>(64, 64, 1, 0);
  auto I1 = ImageLib::create_filled_image<double>(64, 64, 1, 0);

  std::uniform_int_distribution dis(5, 25);

  for(int i = 0; i < nPoints; ++i)
  {
    int px = dis(gen);
    int py = dis(gen);
    I0->set(px, py, 0, 1);
    I1->set(px + delta, py + delta, 0, 1);
  }

  std::string path = "translation/low_density";
  run_and_save_example(I0.get(), I1.get(), path, description);
}

int main(int argc, char** argv)
{
  test1();
  exit(0);
}

/*
def test1():
  description = '''
Translations ought to be exactly achievable even with periodic
boundary conditions. This test verifies that presumption.

It seems that images that are non-smooth on the pixel level causes
divergence. Some binary "smoothing" method may be employed. Instead
of single pixels, small squares or circles could be used.
'''
  nPoints = 30
  delta = 20
  I0 = np.zeros((64,64))
  I1 = I0 + 0
  for i in range(nPoints):
    px = randint(5,25)
    py = randint(5,25)
    I0[px,py] = 1
    I1[px+delta,py+delta] = 1

  subpath = path.join('translation', 'low_density')
  run_and_save_example(I0, I1, subpath, description)
  
  subpath = path.join('translation', 'low_density_smoothed')
  I2 = ndimage.gaussian_filter(I0, sigma=1)
  I3 = ndimage.gaussian_filter(I1, sigma=1)
  run_and_save_example(I2, I3, subpath, description)

  subpath = path.join('translation', 'medium_density')
  nPoints = 200
  for i in range(nPoints):
    px = randint(5,25)
    py = randint(5,25)
    I0[px,py] = 1
    I1[px+delta,py+delta] = 1
  run_and_save_example(I0, I1, subpath, description)

  subpath = path.join('translation', 'medium_density_smoothed')
  I2 = ndimage.gaussian_filter(I0, sigma=1)
  I3 = ndimage.gaussian_filter(I1, sigma=1)
  run_and_save_example(I2, I3, subpath, description)
  
  subpath = path.join('translation', 'high_density')
  nPoints = 400
  for i in range(nPoints):
    px = randint(5,25)
    py = randint(5,25)
    I0[px,py] = 1
    I1[px+delta,py+delta] = 1
  run_and_save_example(I0, I1, subpath, description)
  
  subpath = path.join('translation', 'high_density_smoothed')
  I2 = ndimage.gaussian_filter(I0, sigma=1)
  I3 = ndimage.gaussian_filter(I1, sigma=1)
  run_and_save_example(I2, I3, subpath, description)

  subpath = path.join('translation', 'full_density')
  I0[5:26,5:26] = 1
  I1[(5+delta):(26+delta),(5+delta):(26+delta)] = 1
  run_and_save_example(I0, I1, subpath, description)
  
  subpath = path.join('translation', 'full_density_smoothed')
  I2 = ndimage.gaussian_filter(I0, sigma=1)
  I3 = ndimage.gaussian_filter(I1, sigma=1)
  run_and_save_example(I2, I3, subpath, description)
*/