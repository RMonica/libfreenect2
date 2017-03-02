/*
 * This file is part of the OpenKinect Project. http://www.openkinect.org
 *
 * Copyright (c) 2014 individual OpenKinect contributors. See the CONTRIB file
 * for details.
 *
 * This code is licensed to you under the terms of the Apache License, version
 * 2.0, or, at your option, the terms of the GNU General Public License,
 * version 2.0. See the APACHE20 and GPL2 files for the text of the licenses,
 * or the following URLs:
 * http://www.apache.org/licenses/LICENSE-2.0
 * http://www.gnu.org/licenses/gpl-2.0.txt
 *
 * If you redistribute this file in source form, modified or unmodified, you
 * may:
 *   1) Leave this header intact and distribute it under the same terms,
 *      accompanying it with the APACHE20 and GPL20 files, or
 *   2) Delete the Apache 2.0 clause and accompany it with the GPL2 file, or
 *   3) Delete the GPL v2 clause and accompany it with the APACHE20 file
 * In all cases you must keep the copyright notice intact and include a copy
 * of the CONTRIB file.
 *
 * Binary distributions must follow the binary distribution requirements of
 * either License.
 */

/*******************************************************************************
 * IR SHADOW REMOVAL
 ******************************************************************************/
/*
Author: Riccardo Monica <rmonica[at]ce.unipr.it>
  RIMLab, Department of Information Engineering, University of Parma, Italy
  http://www.rimlab.ce.unipr.it/

  IR Shadow removal filter as explained in the article:
  R. Monica and J. Aleotti. Contour-based next-best view planning from
    point cloud segmentation of unknown objects.
    Autonomous Robots, 2017, DOI:10.1007/s10514-017-9618-0
*/

void kernel prepareShadowFilter(global const float *depths,global uchar * ir_angle_left)
{
  const uint y = get_global_id(0);
  if (y >= 424)
    return;

  const uint base_x = y * 512;

  float current_max_covered_angle = -100.0;

  for (int x = 511; x >= 0; x--)
  {
    const uint i = base_x + x;

    const float depth = depths[i];

    if (depth < MIN_DEPTH)
      continue;

    const float image_offset = -((float)(x) - 256.0);
    const float world_offset = depth * image_offset / FOCAL_LENGTH;
    const float ir_offset = world_offset + IR_EMITTER_DISTANCE;
    const float theta = ir_offset / depth;

    bool ok = true;
    if (current_max_covered_angle < theta)
      current_max_covered_angle = theta;
    else
      ok = false;

    ir_angle_left[i] = ok ? 1 : 0;
  }
}

void kernel applyShadowFilter(global const float *depths,global const uchar * ir_angle_left,global float *filtered)
{
  const uint i = get_global_id(0);

  const uint x = i % 512;
  const uint y = i / 512;

  if (y >= 424)
    return;

  const float depth = depths[i] / 1000.0;

  if (depth < MIN_DEPTH / 1000.0)
  {
    filtered[i] = 0.0f;
    return;
  }

  bool ok = true;
  for (uint di = 0; (di < SHADOW_FILTER_ERODE_LEFT) && (di <= x); di++)
    ok = ok && (ir_angle_left[i - di] != 0);
  for (uint di = 0; (di < SHADOW_FILTER_ERODE_RIGHT) && (x + di < 512); di++)
    ok = ok && (ir_angle_left[i + di] != 0);

  if (x == 0 || y == 0 || x == 511 || y == 423)
  {
    filtered[i] = 0.0f;
    return;
  }

  if (ok && SHADOW_FILTER_MAX_VISUAL_ANGLE < 1.0)
  {
    float centered_x = (float)(x) - 256.0;
    float centered_y = (float)(y) - 212.0;
    float3 w0;
    w0.x = depth * centered_x / FOCAL_LENGTH;
    w0.y = depth * centered_y / FOCAL_LENGTH;
    w0.z = depth;
    #pragma unroll
    for (int dy = -1; dy <= 1; dy++)
      #pragma unroll
      for (int dx = -1; dx <= 1; dx++)
      {
        if (dx && dy)
          continue;
        if (!dx && !dy)
          continue;
        float depth1 = depths[((int)i) + 512 * dy + dx] / 1000.0;
        float3 w1;
        w1.x = depth1 * (centered_x + (float)dx) / FOCAL_LENGTH;
        w1.y = depth1 * (centered_y + (float)dy) / FOCAL_LENGTH;
        w1.z = depth1;
        float3 diff1 = w1 - w0;
        float3 diff2 = w0;
        float cos_angle = fabs(dot(diff1,diff2)) / (length(diff1) * length(diff2));
        if (cos_angle > SHADOW_FILTER_MAX_VISUAL_ANGLE)
          ok = false;
      }
  }

  filtered[i] = ok ? depth * 1000.0 : 0.0f;
}
