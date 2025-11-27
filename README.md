# OpticalCircleFitting
My first idea originated from spiral phase contrast imaging and was inspired by the Hough gradient–based circle fitting method.

During my PhD, my initial plan was to further develop my senior colleague’s work on spiral phase contrast imaging. However, I realized that previous research on spiral phase contrast imaging was already quite comprehensive. So I began to study the mathematical and computer-vision interpretation of spiral phase contrast imaging, and I noticed the phase (i.e., the gradient direction) that had been overlooked during the camera capture process.

This is the first idea that came from my independent thinking, and I’m proud of it. Below I’ll help you get started, and I’ll attach the code and the data from the paper.

1. **First, you need to understand Fourier transform / convolution operations in optical systems.**  
   In optics, this is often called spatial filtering. There are two commonly used spatial-filtering systems.

   - One is implemented by a 4f system made of two lenses.  
   - The other is implemented through superposition of the device’s transmission function with the phase of a single lens.

   In the file “证明4-f和单透镜.md”/"spatial filtering.md", I describe the mechanism in detail and distinguish the differences between these two spatial-filtering systems. The superposition-based system has very high spatial integration, but it has unavoidable defects; the 4f system is more complex spatially, but gives the best spatial-filtering performance.

2. **After understanding optical convolution, you’ll need some computer vision background.**  
   You should know that gradients, as local information, are often obtained by convolving the image with a small kernel. Since optical systems are linear and compatible with convolution, we use the Sobel operator to extract gradients (the Canny operator is more complex and involves nonlinear operations).  

   In computer vision, due to sampling limits, the Sobel operator only extracts gradients in 8 directions: 0°, 45°, 90°, …, 315°. But optical convolution acts in real physical space and is continuous, unlike the discrete operation in CV. Therefore, we extend SobelX + 1i·SobelY to cover all angles over 360°, i.e., a spiral phase, to extract finer gradients. This is also the principle behind spiral phase contrast imaging. The proof is in “菲涅尔衍射-汉克尔变换-螺旋相称成像.md”/"spiarl phase contrasting imaging.md".

From here on is my original thinking.

3. **In (1) I said there are two optical convolution systems, but the superposition system has defects.**  
   It introduces an additional phase that varies from the optical axis center (lens center) to the periphery, and the variation becomes faster farther out. This causes two problems:

   1) Shrinking field of view. Near the optical axis, the extra phase changes slowly—at least within the range of the gradient operator. You can think of our kernel as small: within a small kernel region, the extra term exists but is effectively constant. But far from the axis, even though the kernel is still small, the extra phase grows quadratically and is no longer constant within a local region, so gradient computation fails.  

   2) Distortion of phase / gradient direction. In CV, the intensity of the convolution result is treated as the gradient magnitude, and the phase as the gradient direction. Even near the optical axis where gradient extraction is still valid, that doesn’t mean the extra phase has no effect. It’s only because the kernel is tiny that the extra phase is constant inside one kernel and can be factored out. But different kernels have different extra phases, so the final convolution phase no longer correctly represents the gradient direction.  

   So if we want to use the phase of gradients, we must use a 4f system, to ensure that the spatial filtering / optical convolution phase faithfully represents the gradient direction.

4. **Next is the Hough gradient method.**  
   Its key idea is that the gradient/normal of a circle edge always points to the circle center. So on the gradient map, we use a voter/accumulator: extend the normals (gradients) and count their intersections. The normals from a circular edge intersect at one point, making that point’s weight high.  

   But nonlinear operations are a huge barrier for optical computing. We cannot extend normals only forward without also extending them backward. We also can’t decide whether a normal passes through a point or just grazes past it. This once deeply confused me. I even thought my idea might die here. Damn light: a simple “if-else” in a computer is nearly impossible in optical computing. Optical computing feels like dancing in shackles—basically you can only do Fourier transforms and convolution.

5. **My solution came from thinking about sources and sinks in electrodynamics.**  
   I want the normals on the circle to converge to the circle center. Isn’t that like a sink?  
   Imagine an electromagnetic wave emitted from a source, propagating to a circle/sphere, with its velocity / wave vector along the normal direction from center to periphery. If I reverse the viewpoint: if I give these gradients an extra phase, i.e., rotate the normals by some angle, then they can share the same phase and interfere constructively.  

   This is easy to do: just apply a spiral phase opposite to that in spiral phase contrast imaging to cancel the spiral phase contained in the gradients. This kernel needs to be larger than the circle we want to detect, which feels counter-intuitive. In CV, kernels are usually tiny because we measure local info and large kernels are computationally expensive. But optical convolution’s unique advantage is its computing power!

6. **However, in practice I found I was wrong.**  
   Although I could get a bright spot at the center, it was far too large—almost covering half the circle’s area. I can’t claim I’m doing circle detection and center localization if the result is just “the center is somewhere inside the circle, and this big blob is the center area.” That’s not precise.  
   The reason is that I used a rough convolution to replace the voter in the Hough gradient method. In CV, a normal contributes 1 vote if it passes through the center, and 0 if it only grazes it. But this hard nonlinear decision is impossible in optics. My convolution-based accumulator is not a true voter; when a gradient grazes the center, it still contributes, so the center spot expands.

7. **The inspiration to fix this came from Feynman’s path integral.**  
   Fermat’s principle / least action says light goes from A to B along the path where optical path length is extremal. But Feynman says light is a spherical wave (Huygens principle) and takes all paths, e.g., A → C → B. If we perturb such a path by a tiny δ, the path change is small in geometry but not small compared to the wavelength, so the optical path difference changes dramatically. The vector sum of the perturbed and unperturbed paths becomes like adding random vectors in a chaotic state, yielding nearly zero.  
   This is like air molecules in a box constantly undergoing thermal motion, yet macroscopically we treat air as stationary because velocities cancel. Only paths at the extremum (where optical path difference / phase is stationary) survive and interfere constructively, becoming observable.  
   So if a gradient only grazes the center, can we give it a random chaotic perturbation so it destructively interferes?

8. **Therefore, on top of my convolution-based accumulator, I multiply annulus by annulus outward with 1, −1, 1, −1, …**  
   (equivalently phases 0, π, 0, π, …).  
   For the true center, nothing changes compared to before. But for points around the center, they receive alternating random phases 0/π and enter a chaotic-like state, thus cancelling out by destructive interference.

9. **This also gives an extra surprise.**  
   A ring has an inner circle and outer circle; their gradient directions are opposite (toward vs away from the center). In Hough gradient, that’s fine—we just ignore the outward ones. But in our accumulator, it’s a big issue: when the inner and outer gradients superpose, their opposite phases make them destructively interfere, preventing a bright spot at the center.  
   With the periodic 0/π phase modulation, if the inner and outer gradients fall into two different regions (one in 0, the other in π), then what used to cancel now adds constructively (opposite gradients effectively become aligned).  
   But if both fall into the same 0 region or the same π region, we only get outer minus inner intensity. Since the outer ring is larger, we can still get a bright spot, but weaker.  
   So this relies on luck: instead of guaranteed ring detection, we get about a 50% chance.

10. **From my paper, in practice only two transmitting rings with a π phase difference are enough to achieve circle fitting.**

11. **A future vision: using circle-recognition robots to detect and treat cancer cells.**  
    Current nanorobots either recognize certain target proteins on cancer cell surfaces, or use local pH/temperature cues near tumors to approach and release drugs.  
    The first approach, like targeted drugs, risks selecting for resistant cancer cells and losing effectiveness. The second feels unreliable. Medically, the most convincing way to identify cancer is still biopsy and morphological observation of cell slices—and a key cancer hallmark is an increased nucleus-to-cytoplasm ratio.  
    Some groups have shown that laser irradiation on platinum sheets can bend them, providing propulsion for nanorobots. My idea is to place up-conversion nanoparticles (UCNPs) at the center of the dual-ring circle-fitting device, with 3–4 platinum sheets behind.  
    When irradiated by red or near-infrared light (high tissue penetration), UCNPs emit higher-frequency light as a probe. Because UCNPs are tiny, the emitted light can be treated as spatially coherent. Light reflected from cells passes through the dual-ring circle fitting system and produces a bright spot:

    - the spot position corresponds to the nucleus location;  
    - the spot intensity corresponds to nucleus size (a larger circular nucleus produces more gradients, so more normals accumulate at the center).

    The bright spot hits a specific platinum sheet, bending it and generating propulsion toward the cancerous cell, enabling recognition.

    Of course this sounds very shaky: How to realize spatial filtering in practice? Many cells besides cancer are non-circular. Could the membrane interfere? For internal organ cancers, how to deliver infrared light to the lesion?  
      I’m just tossing out an idea and telling a story. I’m looking forward to researchers’ wisdom.

