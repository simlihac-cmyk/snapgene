"""Generate a simple original PlasmidLab ICO placeholder."""

from __future__ import annotations

import math
import struct
from pathlib import Path


SIZE = 64
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "assets" / "icons" / "plasmidlab.ico"


def main() -> int:
    """Generate the Windows ICO file using only the standard library."""

    pixels = [(0, 0, 0, 0) for _ in range(SIZE * SIZE)]
    _rounded_rect(pixels, 5, 5, 58, 58, 11, (11, 114, 133, 255))
    _circle_stroke(pixels, 32, 32, 17, 5, (231, 245, 255, 255))
    _arc_stroke(pixels, 32, 32, 17, 5, -78, 35, (130, 201, 30, 255))
    _arc_stroke(pixels, 32, 32, 17, 5, 145, 220, (255, 212, 59, 255))
    _filled_circle(pixels, 48, 25, 3, (255, 255, 255, 255))
    _filled_circle(pixels, 20, 20, 3, (255, 255, 255, 255))
    _line(pixels, 25, 32, 39, 32, 4, (255, 255, 255, 255))
    _line(pixels, 32, 25, 32, 39, 4, (255, 255, 255, 255))
    _line(pixels, 21, 49, 43, 49, 3, (231, 245, 255, 255))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_bytes(_ico_bytes(pixels))
    return 0


def _set_pixel(pixels: list[tuple[int, int, int, int]], x: int, y: int, color: tuple[int, int, int, int]) -> None:
    if 0 <= x < SIZE and 0 <= y < SIZE:
        pixels[y * SIZE + x] = color


def _rounded_rect(
    pixels: list[tuple[int, int, int, int]],
    x0: int,
    y0: int,
    x1: int,
    y1: int,
    radius: int,
    color: tuple[int, int, int, int],
) -> None:
    for y in range(y0, y1 + 1):
        for x in range(x0, x1 + 1):
            dx = max(x0 + radius - x, 0, x - (x1 - radius))
            dy = max(y0 + radius - y, 0, y - (y1 - radius))
            if dx * dx + dy * dy <= radius * radius:
                _set_pixel(pixels, x, y, color)


def _filled_circle(
    pixels: list[tuple[int, int, int, int]],
    cx: int,
    cy: int,
    radius: int,
    color: tuple[int, int, int, int],
) -> None:
    for y in range(cy - radius, cy + radius + 1):
        for x in range(cx - radius, cx + radius + 1):
            if (x - cx) ** 2 + (y - cy) ** 2 <= radius * radius:
                _set_pixel(pixels, x, y, color)


def _circle_stroke(
    pixels: list[tuple[int, int, int, int]],
    cx: int,
    cy: int,
    radius: int,
    width: int,
    color: tuple[int, int, int, int],
) -> None:
    for y in range(cy - radius - width, cy + radius + width + 1):
        for x in range(cx - radius - width, cx + radius + width + 1):
            distance = math.hypot(x - cx, y - cy)
            if radius - width / 2 <= distance <= radius + width / 2:
                _set_pixel(pixels, x, y, color)


def _arc_stroke(
    pixels: list[tuple[int, int, int, int]],
    cx: int,
    cy: int,
    radius: int,
    width: int,
    start: float,
    end: float,
    color: tuple[int, int, int, int],
) -> None:
    for y in range(cy - radius - width, cy + radius + width + 1):
        for x in range(cx - radius - width, cx + radius + width + 1):
            distance = math.hypot(x - cx, y - cy)
            if not radius - width / 2 <= distance <= radius + width / 2:
                continue
            angle = math.degrees(math.atan2(y - cy, x - cx)) % 360
            if start % 360 <= angle <= end % 360 or (start > end and (angle >= start or angle <= end)):
                _set_pixel(pixels, x, y, color)


def _line(
    pixels: list[tuple[int, int, int, int]],
    x0: int,
    y0: int,
    x1: int,
    y1: int,
    width: int,
    color: tuple[int, int, int, int],
) -> None:
    steps = max(abs(x1 - x0), abs(y1 - y0), 1)
    for index in range(steps + 1):
        t = index / steps
        x = round(x0 + (x1 - x0) * t)
        y = round(y0 + (y1 - y0) * t)
        _filled_circle(pixels, x, y, width // 2, color)


def _ico_bytes(pixels: list[tuple[int, int, int, int]]) -> bytes:
    header = struct.pack("<HHH", 0, 1, 1)
    bmp_header_size = 40
    xor_size = SIZE * SIZE * 4
    and_mask_size = ((SIZE + 31) // 32 * 4) * SIZE
    image_size = bmp_header_size + xor_size + and_mask_size
    directory = struct.pack("<BBBBHHII", SIZE, SIZE, 0, 0, 1, 32, image_size, 6 + 16)
    bmp_header = struct.pack("<IIIHHIIIIII", 40, SIZE, SIZE * 2, 1, 32, 0, xor_size, 0, 0, 0, 0)
    bgra = bytearray()
    for y in range(SIZE - 1, -1, -1):
        for x in range(SIZE):
            r, g, b, a = pixels[y * SIZE + x]
            bgra.extend((b, g, r, a))
    return header + directory + bmp_header + bytes(bgra) + (b"\x00" * and_mask_size)


if __name__ == "__main__":
    raise SystemExit(main())
