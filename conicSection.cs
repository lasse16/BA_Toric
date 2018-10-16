using UnityEngine;

internal interface conicSection
{
    bool IsInside(Vector3 pointToTest);
    /// <summary>
    /// Draws the conic section in the desired color
    /// </summary>
    /// <param name="color">desired color</param>
    void Draw(Color color);
}