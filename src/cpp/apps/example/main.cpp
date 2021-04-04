#include <filesystem>
#include <iostream>
#include <random>
#include <string>

#include "core/MyArrays.h"

#include "image/Image.h"
#include "image/Image_funcs.h"
#include "image/Image_storage.h"

#include "ddmatch/DiffeoFunctionMatching.h"

using ImageLib::TImage;

namespace fs = std::filesystem;


std::unique_ptr<ImageLib::Image> create_image(const dGrid& grid) {
  int w = grid.cols();
  int h = grid.rows();
  std::unique_ptr<ImageLib::Image> ret = std::make_unique<ImageLib::Image>(w, h, 1);
  uint8_t* dst = ret->data();
  for(int y = 0; y < h; ++y) {
    for(int x = 0; x < w; ++x) {
      int c = static_cast<int>(grid[y][x] * 255.0);
      c = std::min(std::max(c, 0), 255);
      dst[y*w+x] = static_cast<uint8_t>(c & 0xff);
    }
  }
  return ret;
}

bool save_image(const dGrid& grid, const fs::path& filename) {
  auto img = create_image(grid);
  const auto [ok, msg] = ImageLib::save(img.get(), filename.string());
  if (!ok)
    printf("ERROR: %s\n", msg.c_str());
  return ok;
}

dGrid combine_warp(const dGrid& dx, const dGrid& dy) {
  return dx;
}

void run_and_save_example(const dGrid& I0, const dGrid& I1, const std::string& subpath, const std::string& desc)
{
  printf("%s Initializing\n", subpath.c_str());

  bool compute_phi = true;
  std::unique_ptr<DiffeoFunctionMatching> dfm;
  std::string msg;
  double alpha = 0.001;
  double beta  = 0.03;
  double sigma = 0.05;

  //dm = difforma_base.DiffeoFunctionMatching(
  //  source=I0, target=I1,
  //  alpha=0.001, beta=0.03, sigma=0.05
  //)
  std::tie(dfm, msg) = DiffeoFunctionMatching::create(I0, I1, alpha, beta, sigma, compute_phi);
  if (!dfm) {
    printf("ERROR: %s\n", msg.c_str());
    return;
  }

  printf("%s Running\n", subpath.c_str());
  //int num_iters = 1000;
  int num_iters = 1000;
  double epsilon = 0.1;
  dfm->run(num_iters, epsilon);
  
  //if not path.exists(subpath):
  //  os.makedirs(subpath)
  fs::path root_path(subpath);
  printf("%s: Creating plots\n", subpath.c_str());
  fs::create_directories(root_path);

  fs::path overview_path = root_path / "overview";
  fs::create_directories(overview_path);

  //# 2x2 overview plot
  //plt1 = plt.figure(1, figsize=(11.7,9))
  //plt.clf()

  //plt.subplot(2,2,1)
  //plt.imshow(dm.target, cmap='bone', vmin=dm.I0.min(), vmax=dm.I0.max())
  //plt.colorbar()
  //plt.title('Target image')
  save_image(dfm->target(), overview_path/"target.png");

  //plt.subplot(2,2,2)
  //plt.imshow(dm.source, cmap='bone', vmin=dm.I0.min(), vmax=dm.I0.max())
  //plt.colorbar()
  //plt.title('Template image')
  save_image(dfm->source(), overview_path/"template.png");

  //plt.subplot(2,2,3)
  //plt.imshow(dm.I, cmap='bone', vmin=dm.I0.min(), vmax=dm.I0.max())
  //plt.colorbar()
  //plt.title('Warped image')
  save_image(dfm->warped(), overview_path/"warped.png");

  //plt.subplot(2,2,4)
  //# Forward warp    
  //phix = dm.phix
  //phiy = dm.phiy
  //# Uncomment for backward warp
  //#phix = dm.phiinvx
  //#phiy = dm.phiinvy
  //plot_warp(phix, phiy, downsample=4)
  auto warped = combine_warp(dfm->phi_x(), dfm->phi_y());
  save_image(warped, overview_path/"forward_warp.png");
  warped = combine_warp(dfm->phi_inv_x(), dfm->phi_inv_y());
  save_image(warped, overview_path/"backward_warp.png");
  //plt.axis('equal')
  //warplim = [phix.min(), phix.max(), phiy.min(), phiy.max()]
  //warplim[0] = min(warplim[0], warplim[2])
  //warplim[2] = warplim[0]
  //warplim[1] = max(warplim[1], warplim[3])
  //warplim[3] = warplim[1]

  //plt.axis(warplim)
  //plt.gca().invert_yaxis()
  //plt.gca().set_aspect('equal')
  //plt.title('Warp')
  //plt.grid()
  //plt1.savefig(path.join(subpath, 'overview.png'), dpi=300, bbox_inches='tight')

  //# Energy convergence plot
  //plt2 = plt.figure(2, figsize=(8,4.5))
  //plt.clf()
  //plt.plot(dm.E)
  //plt.grid()
  //plt.ylabel('Energy')
  //plt2.savefig(os.path.join(subpath, 'convergence.png'), dpi=150, bbox_inches='tight')

  //# Dedicated warp plot (forward only)
  //plt3 = plt.figure(3, figsize=(10,10))
  //plt.clf()
  //plot_warp(phix, phiy, downsample=4, )
  //plt.axis('equal')
  //warplim = [phix.min(), phix.max(), phiy.min(), phiy.max()]
  //warplim[0] = min(warplim[0], warplim[2])
  //warplim[2] = warplim[0]
  //warplim[1] = max(warplim[1], warplim[3])
  //warplim[3] = warplim[1]

  //plt.axis(warplim)
  //plt.gca().invert_yaxis()
  //plt.gca().set_aspect('equal')
  //plt.title('Warp')
  //plt.axis('off')
  //#plt.grid(color='black')
  //plt3.savefig(path.join(subpath, 'warp.png'), dpi=150, bbox_inches='tight')
  
  //# Output description
  //with open(path.join(subpath, 'description.txt'), 'w') as f:
  //  f.write(subpath)
  //  f.write('\n')
  //  f.write(description)
  //  
  //print('Done at ' + time.asctime() + '\n')
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

void low_density() {
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
  dGrid I0(64,64,0.0);
  dGrid I1(64,64,0.0);

  std::uniform_int_distribution dis(5, 25);

  for(int i = 0; i < nPoints; ++i)
  {
    int px = dis(gen);
    int py = dis(gen);
    I0[py][px] = 1.0;
    I1[py+delta][px+delta] = 1.0;
  }

  std::string path = "translation/low_density";
  run_and_save_example(I0, I1, path, description);
}

void medium_density() {
  //std::random rnd;
  //std::mt19937_64 gen(rnd());
  std::mt19937_64 gen(1234);
  std::string description =
    "Translations ought to be exactly achievable even with periodic\n"
    "boundary conditions. This test verifies that presumption.\n\n"
    "It seems that images that are non-smooth on the pixel level causes divergence.\n"
    "Some binary \"smoothing\" method may be employed. Instead of single pixels,\n"
    "small squares or circles could be used.";
  int nPoints = 200;
  int delta = 20;
  dGrid I0(64, 64, 0.0);
  dGrid I1(64, 64, 0.0);

  std::uniform_int_distribution dis(5, 25);

  for (int i = 0; i < nPoints; ++i)
  {
    int px = dis(gen);
    int py = dis(gen);
    I0[py][px] = 1.0;
    I1[py + delta][px + delta] = 1.0;
  }

  std::string path = "translation/medium_density";
  run_and_save_example(I0, I1, path, description);
}


void high_density() {
  //std::random rnd;
  //std::mt19937_64 gen(rnd());
  std::mt19937_64 gen(1234);
  std::string description =
    "Translations ought to be exactly achievable even with periodic\n"
    "boundary conditions. This test verifies that presumption.\n\n"
    "It seems that images that are non-smooth on the pixel level causes divergence.\n"
    "Some binary \"smoothing\" method may be employed. Instead of single pixels,\n"
    "small squares or circles could be used.";
  int nPoints = 400;
  int delta = 20;
  dGrid I0(64, 64, 0.0);
  dGrid I1(64, 64, 0.0);

  std::uniform_int_distribution dis(5, 25);

  for (int i = 0; i < nPoints; ++i)
  {
    int px = dis(gen);
    int py = dis(gen);
    I0[py][px] = 1.0;
    I1[py + delta][px + delta] = 1.0;
  }

  std::string path = "translation/high_density";
  run_and_save_example(I0, I1, path, description);
}

void full_density() {
  //std::random rnd;
  //std::mt19937_64 gen(rnd());
  std::mt19937_64 gen(1234);
  std::string description =
    "Translations ought to be exactly achievable even with periodic\n"
    "boundary conditions. This test verifies that presumption.\n\n"
    "It seems that images that are non-smooth on the pixel level causes divergence.\n"
    "Some binary \"smoothing\" method may be employed. Instead of single pixels,\n"
    "small squares or circles could be used.";
  int nPoints = 400;
  int delta = 20;
  dGrid I0(64, 64, 0.0);
  dGrid I1(64, 64, 0.0);

  std::uniform_int_distribution dis(5, 25);

  for (int py = 5; py < 25; ++py) {
    for (int px = 5; px < 25; ++px) {
      I0[py][px] = 1.0;
      I1[py + delta][px + delta] = 1.0;
    }
  }
  std::string path = "translation/full_density";
  run_and_save_example(I0, I1, path, description);
}

int main(int argc, char** argv)
{
  low_density();
  medium_density();
  high_density();
  full_density();
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